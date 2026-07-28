#!/bin/bash
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "${REPO_DIR}/scripts/00_util.sh"

work_dir=$1
sample_manifest_tsv=$2
fastqc=$3
MAX_PROCESSOR_USE_PERCENT=$4
MAX_FASTQC_JOBS=${5:-10}
MAX_Q30_STAT_THREADS=${6:-20}
FASTQC_EXTRA_ARGS=${7:-}
multiqc=${8:-multiqc}
seqkit=${9:-seqkit}
conda_bin=${10:-conda}

fastqc_dir="${work_dir}/output/01_fastqc"
mkdir -p "$fastqc_dir"
fastqc_tasks="${fastqc_dir}/fastq_tasks.txt"
: > "$fastqc_tasks"

while IFS= read -r sample_id; do
  [ -z "$sample_id" ] && continue
  read_1=$(find_sample_read "$sample_id" "$sample_manifest_tsv" R1) || exit 1
  read_2=$(find_sample_read "$sample_id" "$sample_manifest_tsv" R2) || exit 1
  if [ -z "$read_1" ] || [ -z "$read_2" ]; then
    echo "Error: raw FASTQ files not found for ${sample_id} using ${sample_manifest_tsv}"
    exit 1
  fi
  printf '%s\n%s\n' "$read_1" "$read_2" >> "$fastqc_tasks"
done < "${work_dir}/sample_list.txt"

num_fastq=$(wc -l < "$fastqc_tasks")
MAX_FASTQC_FILE_JOBS=$(limit_parallel_jobs "$MAX_PROCESSOR_USE_PERCENT" "$((MAX_FASTQC_JOBS * 2))" "$num_fastq")

# Run FastQC on one raw FASTQ. FastQC uses one analysis worker per input file.
run_fastqc() {
  local fastq_file=$1
  local fastq_name output_stem output_html
  local -a extra_args=()
  read_extra_args "$FASTQC_EXTRA_ARGS" extra_args || return 1

  [ -s "$fastq_file" ] || { echo "Error: FASTQ not found or empty: $fastq_file"; return 1; }
  fastq_name=$(basename "$fastq_file")
  output_stem=${fastq_name%.gz}
  output_stem=${output_stem%.fastq}
  output_stem=${output_stem%.fq}
  output_html="${fastqc_dir}/${output_stem}_fastqc.html"
  if [ -s "$output_html" ]; then
    echo -e "\e[37;42mfastqc: ${fastq_name} exists\e[m"
    return 0
  fi

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc: ${fastq_name} processing\e[m"
  "$fastqc" -t 1 "${extra_args[@]}" -o "$fastqc_dir" "$fastq_file"
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc: ${fastq_name} done\e[m"
}

cleanup_run_fastqc() {
  local fastq_name output_stem
  fastq_name=$(basename "$1")
  output_stem=${fastq_name%.gz}
  output_stem=${output_stem%.fastq}
  output_stem=${output_stem%.fq}
  rm -f "${fastqc_dir}/${output_stem}_fastqc.html" "${fastqc_dir}/${output_stem}_fastqc.zip"
}

# Parallel FastQC across individual FASTQ files.
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc: start\e[m"
echo "FastQC resources: ${num_fastq} files, ${MAX_FASTQC_FILE_JOBS} parallel one-thread jobs, CPU budget $(calculate_cpu_budget "$MAX_PROCESSOR_USE_PERCENT")"
parallel_run_sample_list "$fastqc_tasks" "$MAX_FASTQC_FILE_JOBS" run_fastqc fastqc || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc: all done\e[m"

# Aggregate FastQC reports with MultiQC
if [ ! -s "${fastqc_dir}/multiqc_report.html" ]; then
  "$conda_bin" run --no-capture-output -n multiqc "$multiqc" "$fastqc_dir" --outdir "$fastqc_dir"
fi
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc: multiqc done\e[m"

# Seqkit Q30/Q20 stats on raw reads
if [ ! -s "${fastqc_dir}/seqkit_stats.txt" ]; then
  seqkit_inputs=()
  while IFS= read -r sample_id; do
    [ -z "$sample_id" ] && continue
    seqkit_inputs+=("$(find_sample_read "$sample_id" "$sample_manifest_tsv" R1)")
    seqkit_inputs+=("$(find_sample_read "$sample_id" "$sample_manifest_tsv" R2)")
  done < "${work_dir}/sample_list.txt"
  "$conda_bin" run --no-capture-output -n seqkit "$seqkit" stats "${seqkit_inputs[@]}" -aT -j "$MAX_Q30_STAT_THREADS" \
    -o "${fastqc_dir}/seqkit_stats.txt"
fi
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc: seqkit stats done\e[m"
