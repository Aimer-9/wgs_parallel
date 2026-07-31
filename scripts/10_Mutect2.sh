#!/bin/bash
# Stage 10: matched tumor-normal somatic calling with GATK Mutect2.
#
# Positional arguments:
#   1 work_dir; 2 normalized sample manifest; 3 GATK; 4 samtools;
#   5 bcftools; 6 reference build; 7 reference FASTA;
#   8 germline population resource; 9 panel of normals;
#  10 primary_contigs.tsv; 11 maximum concurrent pairs;
#  12 Mutect2 extra args; 13 FilterMutectCalls extra args.
#
# Tumor BAMs are Stage 06 BQSR outputs. A manifest row is eligible only when it
# has normal_bam, an index, a dictionary compatible with all 25 primary contigs,
# and one unambiguous normal sample name matching normal_id when supplied.
# Eligible pairs are recorded in pairs.tsv; validation skips and calling failures
# are recorded separately. Mutect2 is restricted to 1-22, X, Y, and MT. Each
# sample/normal pair is scattered into one Mutect2 task per primary contig;
# chromosome VCFs are gathered before FilterMutectCalls and PASS filtering.
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "${REPO_DIR}/scripts/00_util.sh"

work_dir=$1
sample_manifest_tsv=$2
gatk=$3
samtools=$4
bcftools=$5
hsa_version=$6
ref_fa=$7
germline_resource=$8
panel_of_normals=$9
primary_contigs_tsv=${10}
MAX_MUTECT2_JOBS=${11:-10}
MUTECT2_EXTRA_ARGS=${12:-}
FILTER_MUTECT_EXTRA_ARGS=${13:-}

mutect2_dir="${work_dir}/output/10_Mutect2"
bqsr_dir="${work_dir}/output/06_bqsr"
pair_manifest="${mutect2_dir}/pairs.tsv"
pending_pair_manifest="${mutect2_dir}/pending_pairs.tsv"
interval_manifest="${mutect2_dir}/intervals.tsv"
skipped_pairs="${mutect2_dir}/skipped_pairs.tsv"
failed_pairs="${mutect2_dir}/failed_pairs.tsv"
success_samples="${mutect2_dir}/success_samples.txt"
mkdir -p "$mutect2_dir"
printf 'sample_id\tnormal_bam\treason\n' > "$skipped_pairs"
: > "$pair_manifest"
: > "$pending_pair_manifest"
: > "$interval_manifest"
: > "$success_samples"

# Record non-runnable rows without stopping valid tumor-normal pairs.
skip_pair() {
  printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "$skipped_pairs"
  echo "Mutect2 skip: $1 - $3"
}

# Confirm that every selected reference contig exists in the normal BAM header
# with the same length. This catches hg19/hg38 and chr/non-chr mismatches early.
normal_has_primary_dictionary() {
  local normal_bam=$1
  local header
  header=$("$samtools" view -H "$normal_bam" 2>/dev/null) || return 1
  while IFS=$'\t' read -r order logical contig length; do
    if ! awk -v want="$contig" -v want_len="$length" -F '\t' '
      $1=="@SQ" {
        sn=""; ln="";
        for (i=2; i<=NF; i++) {
          if ($i ~ /^SN:/) sn=substr($i,4)
          if ($i ~ /^LN:/) ln=substr($i,4)
        }
        if (sn==want && ln==want_len) found=1
      }
      END { exit !found }
    ' <<< "$header"; then
      return 1
    fi
  done < "$primary_contigs_tsv"
}

# Validate normal files and sample identities before launching any caller jobs.
while IFS= read -r manifest_row; do
  normalized_row=${manifest_row//$'\t'/$'\x1f'}
  IFS=$'\x1f' read -r sample_id sample_prefix sample_ref raw_dir normal_id normal_bam <<< "$normalized_row"
  [ -z "$sample_id" ] || [ -z "$normal_bam" ] && continue
  if [ ! -s "$normal_bam" ]; then
    skip_pair "$sample_id" "$normal_bam" "normal BAM not found or empty"
    continue
  fi
  if [ ! -s "${normal_bam}.bai" ] && [ ! -s "${normal_bam%.bam}.bai" ] && [ ! -s "${normal_bam}.csi" ]; then
    skip_pair "$sample_id" "$normal_bam" "normal BAM index not found"
    continue
  fi
  tumor_bam="${bqsr_dir}/${sample_id}.sorted.markdup.bqsr.bam"
  if [ ! -s "$tumor_bam" ]; then
    skip_pair "$sample_id" "$normal_bam" "tumor BQSR BAM not found"
    continue
  fi
  if ! normal_has_primary_dictionary "$normal_bam"; then
    skip_pair "$sample_id" "$normal_bam" "normal BAM dictionary is incompatible with the selected reference"
    continue
  fi
  sample_name_file="${mutect2_dir}/.${sample_id}.normal_sample.txt"
  if ! "$gatk" GetSampleName -I "$normal_bam" -O "$sample_name_file" >/dev/null 2>&1; then
    skip_pair "$sample_id" "$normal_bam" "cannot determine normal BAM sample name"
    rm -f "$sample_name_file"
    continue
  fi
  normal_sample=$(tr -d '\r\n' < "$sample_name_file")
  rm -f "$sample_name_file"
  if [ -z "$normal_sample" ] || [[ "$normal_sample" == *$'\t'* ]]; then
    skip_pair "$sample_id" "$normal_bam" "normal BAM must contain one unambiguous sample name"
    continue
  fi
  if [ -n "$normal_id" ] && [ "$normal_id" != "$normal_sample" ]; then
    skip_pair "$sample_id" "$normal_bam" "normal_id '${normal_id}' does not match BAM sample '${normal_sample}'"
    continue
  fi
  if [ "$normal_sample" = "$sample_id" ]; then
    skip_pair "$sample_id" "$normal_bam" "tumor and normal sample names are identical"
    continue
  fi
  printf '%s\t%s\t%s\n' "$sample_id" "$normal_bam" "$normal_sample" >> "$pair_manifest"
  pass_vcf="${mutect2_dir}/${sample_id}.mutect2.filtered.PASS.vcf"
  if [ -s "$pass_vcf" ]; then
    echo "Mutect2: ${sample_id} final PASS VCF exists; skipping calling"
    while IFS=$'\t' read -r order logical contig length; do
      rm -f "${mutect2_dir}/${sample_id}.mutect2.${contig}.unfiltered.vcf"*
    done < "$primary_contigs_tsv"
    continue
  fi
  printf '%s\t%s\t%s\n' "$sample_id" "$normal_bam" "$normal_sample" >> "$pending_pair_manifest"
  while IFS=$'\t' read -r order logical contig length; do
    printf '%s\t%s\t%s\t%s\n' "$sample_id" "$normal_bam" "$normal_sample" "$contig" >> "$interval_manifest"
  done < "$primary_contigs_tsv"
done < "$sample_manifest_tsv"

run_mutect2_interval() {
  local task=$1
  local sample_id normal_bam normal_sample contig
  IFS=$'\t' read -r sample_id normal_bam normal_sample contig <<< "$task"
  local tumor_bam="${bqsr_dir}/${sample_id}.sorted.markdup.bqsr.bam"
  local interval_vcf="${mutect2_dir}/${sample_id}.mutect2.${contig}.unfiltered.vcf"
  local -a mutect_extra=()
  read_extra_args "$MUTECT2_EXTRA_ARGS" mutect_extra || return 1

  if [ ! -s "$interval_vcf" ]; then
    echo "Mutect2: ${sample_id} ${contig} processing"
    "$gatk" Mutect2 \
      --verbosity ERROR \
      -R "$ref_fa" \
      -I "$tumor_bam" \
      -I "$normal_bam" \
      -normal "$normal_sample" \
      --germline-resource "$germline_resource" \
      --panel-of-normals "$panel_of_normals" \
      -L "$contig" \
      "${mutect_extra[@]}" \
      -O "$interval_vcf" || return 1
  fi
  [ -s "$interval_vcf" ]
}

cleanup_run_mutect2_interval() {
  local sample_id normal_bam normal_sample contig
  IFS=$'\t' read -r sample_id normal_bam normal_sample contig <<< "$1"
  rm -f "${mutect2_dir}/${sample_id}.mutect2.${contig}.unfiltered.vcf"*
}

finalize_mutect2_pair() {
  local task=$1
  local sample_id normal_bam normal_sample
  IFS=$'\t' read -r sample_id normal_bam normal_sample <<< "$task"
  local unfiltered_vcf="${mutect2_dir}/${sample_id}.mutect2.unfiltered.vcf"
  local merged_stats="${unfiltered_vcf}.stats"
  local filtered_vcf="${mutect2_dir}/${sample_id}.mutect2.filtered.vcf"
  local pass_vcf="${mutect2_dir}/${sample_id}.mutect2.filtered.PASS.vcf"
  local -a filter_extra=() gather_inputs=() stats_inputs=() interval_outputs=()
  read_extra_args "$FILTER_MUTECT_EXTRA_ARGS" filter_extra || return 1

  while IFS=$'\t' read -r order logical contig length; do
    interval_vcf="${mutect2_dir}/${sample_id}.mutect2.${contig}.unfiltered.vcf"
    [ -s "$interval_vcf" ] || return 1
    [ -s "${interval_vcf}.stats" ] || return 1
    gather_inputs+=("-I" "$interval_vcf")
    stats_inputs+=("-stats" "${interval_vcf}.stats")
    interval_outputs+=("$interval_vcf" "${interval_vcf}.idx" "${interval_vcf}.stats")
  done < "$primary_contigs_tsv"

  if [ ! -s "$unfiltered_vcf" ]; then
    echo "GatherVcfs: ${sample_id} processing"
    "$gatk" GatherVcfs "${gather_inputs[@]}" -O "$unfiltered_vcf" || return 1
  fi
  if [ ! -s "$merged_stats" ]; then
    echo "MergeMutectStats: ${sample_id} processing"
    "$gatk" MergeMutectStats "${stats_inputs[@]}" -O "$merged_stats" || return 1
  fi
  if [ -s "$unfiltered_vcf" ] && [ -s "$merged_stats" ]; then
    rm -f "${interval_outputs[@]}"
  fi
  if [ ! -s "$filtered_vcf" ]; then
    echo "FilterMutectCalls: ${sample_id} processing"
    "$gatk" FilterMutectCalls \
      --verbosity ERROR \
      -R "$ref_fa" \
      -V "$unfiltered_vcf" \
      --stats "$merged_stats" \
      "${filter_extra[@]}" \
      -O "$filtered_vcf" || return 1
  fi
  if [ ! -s "$pass_vcf" ]; then
    "$bcftools" view -f PASS "$filtered_vcf" -o "$pass_vcf" || return 1
  fi
  [ -s "$pass_vcf" ]
}

cleanup_finalize_mutect2_pair() {
  local sample_id normal_bam normal_sample
  IFS=$'\t' read -r sample_id normal_bam normal_sample <<< "$1"
  rm -f "${mutect2_dir}/${sample_id}.mutect2.unfiltered.vcf"* \
    "${mutect2_dir}/${sample_id}.mutect2.filtered.vcf"* \
    "${mutect2_dir}/${sample_id}.mutect2.filtered.PASS.vcf"*
}

mutect_status=0
if [ -s "$pair_manifest" ]; then
  echo "Mutect2: start"
  if [ -s "$interval_manifest" ]; then
    parallel_run_sample_list "$interval_manifest" "$MAX_MUTECT2_JOBS" run_mutect2_interval mutect2_interval never || mutect_status=$?
  fi
  if [ -s "$pending_pair_manifest" ]; then
    parallel_run_sample_list "$pending_pair_manifest" "$MAX_MUTECT2_JOBS" finalize_mutect2_pair mutect2_finalize never || mutect_status=$?
  fi
else
  echo "Mutect2: no valid tumor-normal pairs"
fi

# Reconcile expected pair outputs after Parallel finishes. Samples with a final
# gathered and filtered VCF are passed to Stage 11 somatic ANNOVAR annotation.
printf 'sample_id\tnormal_bam\treason\n' > "$failed_pairs"
while IFS=$'\t' read -r sample_id normal_bam normal_sample; do
  [ -z "$sample_id" ] && continue
  if [ -s "${mutect2_dir}/${sample_id}.mutect2.filtered.vcf" ]; then
    echo "$sample_id" >> "$success_samples"
  else
    printf '%s\t%s\tcalling or filtering failed\n' "$sample_id" "$normal_bam" >> "$failed_pairs"
    mutect_status=1
  fi
done < "$pair_manifest"

echo "Mutect2: all runnable pairs attempted"
exit "$mutect_status"
