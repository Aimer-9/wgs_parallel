#!/bin/bash
# Stage 07: chromosome-scattered GATK HaplotypeCaller in GVCF mode.
#
# Positional arguments:
#   1 work_dir; 2 reference build; 3 GATK; 4 reference FASTA;
#   5 maximum concurrent chromosome tasks; 6 primary_contigs.tsv;
#   7 native PairHMM threads/task; 8 keep scatter GVCFs (yes/no);
#   9 additional HaplotypeCaller arguments.
#
# Each incomplete sample is expanded into exactly 25 tasks: chromosomes 1-22,
# X, Y, and mitochondrial DNA using the reference's actual contig spelling.
# Shards are written beneath output/07_HaplotypeCaller/scatter/<sample>, then
# gathered in deterministic contig order into <sample>.g.vcf.gz. Samples with a
# complete indexed final GVCF bypass scatter on restart. Interrupted checkpoints
# remove only the affected shard or final gather output before retrying.
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "${REPO_DIR}/scripts/00_util.sh"

work_dir=$1
hsa_version=$2
gatk=$3
ref_fa=$4
MAX_HC_JOBS=${5:-10}
primary_contigs_tsv=$6
HC_NATIVE_THREADS=${7:-1}
HC_KEEP_SCATTER=${8:-no}
HC_EXTRA_ARGS=${9:-}

bqsr_dir="${work_dir}/output/06_bqsr"
haplotype_dir="${work_dir}/output/07_HaplotypeCaller"
scatter_root="${haplotype_dir}/scatter"
task_file="${scatter_root}/tasks.tsv"
mkdir -p "$scatter_root"
: > "$task_file"

# Do not recreate disposable shards when the final indexed GVCF already proves
# that a prior scatter/gather completed successfully.
while IFS= read -r sample_id; do
  [ -z "$sample_id" ] && continue
  final_gvcf="${haplotype_dir}/${sample_id}.g.vcf.gz"
  if [ -s "$final_gvcf" ] && [ -s "${final_gvcf}.tbi" ]; then
    continue
  fi
  while IFS=$'\t' read -r order logical contig length; do
    printf '%s\t%s\t%s\t%s\n' "$sample_id" "$order" "$logical" "$contig" >> "$task_file"
  done < "$primary_contigs_tsv"
done < "${work_dir}/sample_list.txt"

run_haplotype_shard() {
  local task=$1
  local sample_id order logical contig
  IFS=$'\t' read -r sample_id order logical contig <<< "$task"
  local bqsr_bam="${bqsr_dir}/${sample_id}.sorted.markdup.bqsr.bam"
  local sample_dir="${scatter_root}/${sample_id}"
  local gvcf="${sample_dir}/${order}.${logical}.g.vcf.gz"
  local -a extra_args=()
  read_extra_args "$HC_EXTRA_ARGS" extra_args || return 1

  if [ ! -s "$bqsr_bam" ] || { [ ! -s "${bqsr_bam}.bai" ] && [ ! -s "${bqsr_bam%.bam}.bai" ]; }; then
    echo "Error: indexed BQSR BAM not found for ${sample_id}"
    return 1
  fi
  mkdir -p "$sample_dir"
  if [ -s "$gvcf" ] && [ -s "${gvcf}.tbi" ]; then
    echo "HaplotypeCaller: ${sample_id} ${contig} exists"
    return 0
  fi
  # A complete shard without an index can be repaired without recalling it.
  if [ -s "$gvcf" ]; then
    "$gatk" IndexFeatureFile -I "$gvcf" || return 1
    [ -s "${gvcf}.tbi" ] && return 0
  fi

  date +"%Y-%m-%d %H:%M:%S"
  echo "HaplotypeCaller: ${sample_id} ${contig} processing"
  "$gatk" HaplotypeCaller \
    --verbosity WARNING \
    -R "$ref_fa" \
    -I "$bqsr_bam" \
    -L "$contig" \
    --native-pair-hmm-threads "$HC_NATIVE_THREADS" \
    -ERC GVCF \
    "${extra_args[@]}" \
    -O "$gvcf" || return 1
  [ -s "$gvcf" ] && [ -s "${gvcf}.tbi" ]
}

cleanup_run_haplotype_shard() {
  local sample_id order logical contig
  IFS=$'\t' read -r sample_id order logical contig <<< "$1"
  local gvcf="${scatter_root}/${sample_id}/${order}.${logical}.g.vcf.gz"
  rm -f "$gvcf" "${gvcf}.tbi" "${gvcf}.idx"
}

gather_haplotype_sample() {
  local sample_id=$1
  local final_gvcf="${haplotype_dir}/${sample_id}.g.vcf.gz"
  local sample_dir="${scatter_root}/${sample_id}"
  local order logical contig length shard
  local -a inputs=()

  if [ -s "$final_gvcf" ] && [ -s "${final_gvcf}.tbi" ]; then
    echo "HaplotypeCaller gather: ${sample_id} exists"
    return 0
  fi
  # GatherVcfs receives shards in primary_contigs.tsv order, not filesystem or
  # lexical order, so chr10 cannot precede chr2 accidentally.
  while IFS=$'\t' read -r order logical contig length; do
    shard="${sample_dir}/${order}.${logical}.g.vcf.gz"
    if [ ! -s "$shard" ] || [ ! -s "${shard}.tbi" ]; then
      echo "Error: HaplotypeCaller shard missing for ${sample_id}: ${logical}"
      return 1
    fi
    inputs+=("-I" "$shard")
  done < "$primary_contigs_tsv"

  echo "HaplotypeCaller gather: ${sample_id} processing"
  "$gatk" GatherVcfs "${inputs[@]}" -O "$final_gvcf" || return 1
  if [ ! -s "${final_gvcf}.tbi" ]; then
    "$gatk" IndexFeatureFile -I "$final_gvcf" || return 1
  fi
  # Shards are restart artifacts after a successful gather and may be removed
  # to save disk unless the user explicitly requests retention.
  if [ "$HC_KEEP_SCATTER" = "no" ]; then
    rm -f "${sample_dir}"/*.g.vcf.gz "${sample_dir}"/*.g.vcf.gz.tbi
    rmdir "$sample_dir" 2>/dev/null || true
  fi
}

cleanup_gather_haplotype_sample() {
  local final_gvcf="${haplotype_dir}/$1.g.vcf.gz"
  rm -f "$final_gvcf" "${final_gvcf}.tbi" "${final_gvcf}.idx"
}

date +"%Y-%m-%d %H:%M:%S"
echo "HaplotypeCaller scatter: start"
if [ -s "$task_file" ]; then
  parallel_run_sample_list "$task_file" "$MAX_HC_JOBS" run_haplotype_shard haplotype_scatter || exit 1
else
  echo "HaplotypeCaller scatter: all final gVCFs already complete"
fi
echo "HaplotypeCaller gather: start"
parallel_run_sample_list "${work_dir}/sample_list.txt" "$MAX_HC_JOBS" gather_haplotype_sample haplotype_gather || exit 1
date +"%Y-%m-%d %H:%M:%S"
echo "HaplotypeCaller: all done"
