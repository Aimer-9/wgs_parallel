#!/bin/bash
# Stage 06: GATK base quality score recalibration (BQSR).
#
# Positional arguments:
#   1 work_dir; 2 reference build; 3 GATK; 4 samtools; 5 reference FASTA;
#   6 known indels; 7 dbSNP; 8 Mills indels; 9 CPU fraction;
#  10 maximum concurrent samples; 11 BaseRecalibrator extra arguments.
#
# Input BAMs come from output/05_markduplicated. For every sample this stage
# builds a recalibration report, applies it to a new BAM, and indexes that BAM
# under output/06_bqsr. Checkpoint recovery conservatively removes the report,
# BAM, and indexes for an interrupted sample so dependent work is regenerated.
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "${REPO_DIR}/scripts/00_util.sh"

work_dir=$1
hsa_version=$2
gatk=$3
samtools=$4
ref_fa=$5
known_indels=$6
dbsnp=$7
Mills=$8
MAX_PROCESSOR_USE_PERCENT=$9
MAX_BQSR_JOBS=${10:-10}
BQSR_EXTRA_ARGS=${11:-}

num_samples=$(wc -l < "${work_dir}/sample_list.txt")
nproc=$(nproc)
[ "$MAX_BQSR_JOBS" -gt "$num_samples" ] && MAX_BQSR_JOBS=$num_samples

threads_per_sample=$(calculate_threads "$MAX_PROCESSOR_USE_PERCENT" "$MAX_BQSR_JOBS")

markdup_dir="${work_dir}/output/05_markduplicated"
bqsr_dir="${work_dir}/output/06_bqsr"
mkdir -p "$bqsr_dir"

# BaseRecalibrator, ApplyBQSR, and index for one sample. ApplyBQSR is rerun when
# an existing BAM is smaller than the configured minimum output size.
run_bqsr() {
  local sample_id=$1
  local markdup_bam="${markdup_dir}/${sample_id}.sorted.markdup.bam"
  local bqsr_report="${bqsr_dir}/${sample_id}.sorted.markdup.bqsr.report"
  local bqsr_bam="${bqsr_dir}/${sample_id}.sorted.markdup.bqsr.bam"
  local -a extra_args=()
  read_extra_args "$BQSR_EXTRA_ARGS" extra_args || return 1

  if [ ! -s "$markdup_bam" ]; then
    echo "Error: ${markdup_bam} not found — markdup step may have failed"; return 1
  fi

  # BaseRecalibrator
  if [ ! -s "$bqsr_report" ]; then
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mBQSR recal: ${sample_id} processing\e[m"
    "$gatk" BaseRecalibrator \
      --verbosity WARNING \
      -I "$markdup_bam" \
      -R "$ref_fa" \
      -O "$bqsr_report" \
      --known-sites "$dbsnp" \
      --known-sites "$known_indels" \
      --known-sites "$Mills" "${extra_args[@]}"
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mBQSR recal: ${sample_id} done\e[m"
  else
    echo -e "\e[37;42mBQSR recal: ${sample_id} exists\e[m"
  fi

  if [ ! -s "$bqsr_report" ]; then
    echo "Error: BQSR report missing for ${sample_id} — BaseRecalibrator may have failed"; return 1
  fi

  if [ ! -s "$bqsr_bam" ]; then
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mApplyBQSR: ${sample_id} processing\e[m"
    "$gatk" ApplyBQSR \
      --verbosity WARNING \
      -I "$markdup_bam" \
      -R "$ref_fa" \
      --bqsr-recal-file "$bqsr_report" \
      -O "$bqsr_bam"
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mApplyBQSR: ${sample_id} done\e[m"
  else
    echo -e "\e[37;42mApplyBQSR: ${sample_id} exists\e[m"
  fi

  if [ ! -s "$bqsr_bam" ]; then
    echo "Error: BQSR BAM missing for ${sample_id} — ApplyBQSR may have failed"; return 1
  fi

  # Index
  if [ ! -s "${bqsr_bam}.bai" ]; then
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mBQSR index: ${sample_id} processing\e[m"
    "$samtools" index -@ "$threads_per_sample" "$bqsr_bam"
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mBQSR index: ${sample_id} done\e[m"
  else
    echo -e "\e[37;42mBQSR index: ${sample_id} exists\e[m"
  fi
}

cleanup_run_bqsr() {
  local sample_id=$1
  local bqsr_bam="${bqsr_dir}/${sample_id}.sorted.markdup.bqsr.bam"
  rm -f "${bqsr_dir}/${sample_id}.sorted.markdup.bqsr.report" \
    "$bqsr_bam" "${bqsr_bam}.bai" "${bqsr_bam%.bam}.bai"
}

date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mBQSR: start\e[m"
parallel_run_sample_list "${work_dir}/sample_list.txt" "$MAX_BQSR_JOBS" run_bqsr bqsr || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mBQSR: all done\e[m"
