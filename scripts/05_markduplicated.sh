#!/bin/bash
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "${REPO_DIR}/scripts/00_util.sh"

work_dir=$1
hsa_version=$2
gatk=$3
samtools=$4
MAX_PROCESSOR_USE_PERCENT=$5
MAX_MARKDUP_JOBS=${6:-10}
MARKDUP_EXTRA_ARGS=${7:-}

num_samples=$(wc -l < "${work_dir}/sample_list.txt")
nproc=$(nproc)
[ "$MAX_MARKDUP_JOBS" -gt "$num_samples" ] && MAX_MARKDUP_JOBS=$num_samples

threads_per_sample=$(calculate_threads "$MAX_PROCESSOR_USE_PERCENT" "$MAX_MARKDUP_JOBS")

bwa_dir="${work_dir}/output/04_bwa"
markdup_dir="${work_dir}/output/05_markduplicated"
if [ ! -d "$markdup_dir" ]; then mkdir -p "$markdup_dir"; fi

# Mark duplicates and index for one sample.
run_markdup() {
  local sample_id=$1
  local sorted_bam="${bwa_dir}/${sample_id}.sorted.bam"
  local markdup_bam="${markdup_dir}/${sample_id}.sorted.markdup.bam"
  local metrics="${markdup_dir}/${sample_id}.marked_dup_metrics.txt"
  local -a extra_args=()
  read_extra_args "$MARKDUP_EXTRA_ARGS" extra_args || return 1

  if [ ! -s "$sorted_bam" ]; then
    echo "Error: ${sorted_bam} not found — bwa/sort step may have failed"; return 1
  fi

  if [ ! -s "$markdup_bam" ]; then
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mmarkdup: ${sample_id} processing\e[m"
    "$gatk" MarkDuplicates \
      --VERBOSITY WARNING \
      -I "$sorted_bam" \
      -O "$markdup_bam" \
      -M "$metrics" "${extra_args[@]}"
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mmarkdup: ${sample_id} done\e[m"
  else
    echo -e "\e[37;42mmarkdup: ${sample_id} exists\e[m"
  fi

  if [ ! -s "$markdup_bam" ]; then
    echo "Error: markdup output missing for ${sample_id} — MarkDuplicates may have failed"; return 1
  fi

  if [ ! -s "${markdup_bam}.bai" ]; then
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mmarkdup index: ${sample_id} processing\e[m"
    "$samtools" index -@ "$threads_per_sample" "$markdup_bam"
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mmarkdup index: ${sample_id} done\e[m"
  else
    echo -e "\e[37;42mmarkdup index: ${sample_id} exists\e[m"
  fi
}

cleanup_run_markdup() {
  local sample_id=$1
  local markdup_bam="${markdup_dir}/${sample_id}.sorted.markdup.bam"
  rm -f "$markdup_bam" "${markdup_bam}.bai" \
    "${markdup_bam%.bam}.bai" "${markdup_dir}/${sample_id}.marked_dup_metrics.txt"
}

date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mmarkdup: start\e[m"
parallel_run_sample_list "${work_dir}/sample_list.txt" "$MAX_MARKDUP_JOBS" run_markdup markdup || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mmarkdup: all done\e[m"
