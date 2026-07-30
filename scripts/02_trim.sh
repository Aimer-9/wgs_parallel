#!/bin/bash
# Stage 02: paired-end adapter and quality trimming with Trimmomatic.
#
# Positional arguments:
#   1 work_dir; 2 sample_manifest_tsv; 3 Trimmomatic installation directory;
#   4 CPU fraction; 5 maximum concurrent samples; 6 LEADING quality;
#   7 TRAILING quality; 8 MINLEN; 9 SLIDINGWINDOW; 10 extra arguments;
#  11 Java executable.
#
# For each sample, raw R1/R2 files are resolved from its manifest raw_dir and
# four gzip files are created in output/02_trim: paired and unpaired output for
# each mate. Downstream alignment consumes only the paired files. Trimmomatic
# threads are divided across concurrent sample jobs. On interrupted recovery,
# all four files for that sample are removed before the task is retried.
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "${REPO_DIR}/scripts/00_util.sh"

work_dir=$1
sample_manifest_tsv=$2
Trimmomatic_path=$3
MAX_PROCESSOR_USE_PERCENT=$4
MAX_TRIM_JOBS=${5:-10}
TRIM_LEADING=${6:-10}
TRIM_TRAILING=${7:-10}
TRIM_MINLEN=${8:-25}
TRIM_SLIDINGWINDOW=${9:-4:15}
TRIM_EXTRA_ARGS=${10:-}
java=${11:-java}

num_samples=$(wc -l < "${work_dir}/sample_list.txt")
nproc=$(nproc)
if [ "$MAX_TRIM_JOBS" -gt "$num_samples" ]; then
  MAX_TRIM_JOBS=$num_samples
fi
threads_per_sample=$(calculate_threads "$MAX_PROCESSOR_USE_PERCENT" "$MAX_TRIM_JOBS")

trim_dir="${work_dir}/output/02_trim"
mkdir -p "$trim_dir"

# Trim one sample with Trimmomatic PE mode.
# Skips if paired output already exists. Fails if input reads are not found.
run_trim() {
  local sample_id=$1
  local read_1 read_2
  local -a extra_args=()
  read_extra_args "$TRIM_EXTRA_ARGS" extra_args || return 1

  [ -n "$sample_id" ] || { echo "Error: empty sample_id passed to trimming"; return 1; }
  read_1=$(find_sample_read "$sample_id" "$sample_manifest_tsv" R1) || return 1
  read_2=$(find_sample_read "$sample_id" "$sample_manifest_tsv" R2) || return 1

  if [ -z "$read_1" ] || [ -z "$read_2" ]; then
    echo "Error: raw FASTQ files not found for ${sample_id} using ${sample_manifest_tsv}"
    return 1
  fi

  if [ -s "${trim_dir}/${sample_id}_1_paired.fq.gz" ]; then
    echo -e "\e[37;42mtrim: ${sample_id} exists\e[m"
    return 0
  fi

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mtrim: ${sample_id} processing\e[m"
  "$java" -jar "${Trimmomatic_path}/trimmomatic-0.39.jar" PE \
    "$read_1" "$read_2" \
    "${trim_dir}/${sample_id}_1_paired.fq.gz" \
    "${trim_dir}/${sample_id}_1_unpaired.fq.gz" \
    "${trim_dir}/${sample_id}_2_paired.fq.gz" \
    "${trim_dir}/${sample_id}_2_unpaired.fq.gz" \
    ILLUMINACLIP:"${Trimmomatic_path}/adapters/TruSeq3-PE.fa":2:30:10 \
    "LEADING:${TRIM_LEADING}" "TRAILING:${TRIM_TRAILING}" "MINLEN:${TRIM_MINLEN}" "SLIDINGWINDOW:${TRIM_SLIDINGWINDOW}" \
    -phred33 \
    -threads "$threads_per_sample" "${extra_args[@]}"
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mtrim: ${sample_id} done\e[m"
}

cleanup_run_trim() {
  local sample_id=$1
  rm -f "${trim_dir}/${sample_id}_1_paired.fq.gz" \
    "${trim_dir}/${sample_id}_1_unpaired.fq.gz" \
    "${trim_dir}/${sample_id}_2_paired.fq.gz" \
    "${trim_dir}/${sample_id}_2_unpaired.fq.gz"
}

# Parallel trimming across all samples
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mtrim: start\e[m"
parallel_run_sample_list "${work_dir}/sample_list.txt" "$MAX_TRIM_JOBS" run_trim trim || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mtrim: all done\e[m"
