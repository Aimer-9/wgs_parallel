#!/bin/bash
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "${REPO_DIR}/scripts/00_util.sh"

work_dir=$1
fastqc=$2
MAX_PROCESSOR_USE_PERCENT=$3
MAX_FASTQC_JOBS=${4:-10}
MAX_Q30_STAT_THREADS=${5:-20}
FASTQC_EXTRA_ARGS=${6:-}
multiqc=${7:-multiqc}
seqkit=${8:-seqkit}
conda_bin=${9:-conda}

trim_dir="${work_dir}/output/02_trim"
fastqc_trim_dir="${work_dir}/output/03_fastqc_trim"
mkdir -p "$fastqc_trim_dir"
fastqc_tasks="${fastqc_trim_dir}/fastq_tasks.txt"
: > "$fastqc_tasks"

while IFS= read -r sample_id; do
  [ -z "$sample_id" ] && continue
  paired_1="${trim_dir}/${sample_id}_1_paired.fq.gz"
  paired_2="${trim_dir}/${sample_id}_2_paired.fq.gz"
  if [ ! -s "$paired_1" ] || [ ! -s "$paired_2" ]; then
    echo "Error: trimmed FASTQ not found for ${sample_id} - run 02_trim.sh first"
    exit 1
  fi
  printf '%s\n%s\n' "$paired_1" "$paired_2" >> "$fastqc_tasks"
done < "${work_dir}/sample_list.txt"

num_fastq=$(wc -l < "$fastqc_tasks")
MAX_FASTQC_FILE_JOBS=$(limit_parallel_jobs "$MAX_PROCESSOR_USE_PERCENT" "$((MAX_FASTQC_JOBS * 2))" "$num_fastq")

# Run FastQC on one trimmed FASTQ.
run_fastqc_trim() {
  local fastq_file=$1
  local fastq_name output_stem output_html
  local -a extra_args=()
  read_extra_args "$FASTQC_EXTRA_ARGS" extra_args || return 1

  [ -s "$fastq_file" ] || { echo "Error: FASTQ not found or empty: $fastq_file"; return 1; }
  fastq_name=$(basename "$fastq_file")
  output_stem=${fastq_name%.gz}
  output_stem=${output_stem%.fastq}
  output_stem=${output_stem%.fq}
  output_html="${fastqc_trim_dir}/${output_stem}_fastqc.html"
  if [ -s "$output_html" ]; then
    echo -e "\e[37;42mfastqc_trim: ${fastq_name} exists\e[m"
    return 0
  fi

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc_trim: ${fastq_name} processing\e[m"
  "$fastqc" -t 1 "${extra_args[@]}" -o "$fastqc_trim_dir" "$fastq_file"
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc_trim: ${fastq_name} done\e[m"
}

cleanup_run_fastqc_trim() {
  local fastq_name output_stem
  fastq_name=$(basename "$1")
  output_stem=${fastq_name%.gz}
  output_stem=${output_stem%.fastq}
  output_stem=${output_stem%.fq}
  rm -f "${fastqc_trim_dir}/${output_stem}_fastqc.html" \
    "${fastqc_trim_dir}/${output_stem}_fastqc.zip"
}

# Parallel FastQC across individual FASTQ files.
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc_trim: start\e[m"
echo "FastQC resources: ${num_fastq} files, ${MAX_FASTQC_FILE_JOBS} parallel one-thread jobs, CPU budget $(calculate_cpu_budget "$MAX_PROCESSOR_USE_PERCENT")"
parallel_run_sample_list "$fastqc_tasks" "$MAX_FASTQC_FILE_JOBS" run_fastqc_trim fastqc_trim || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc_trim: all done\e[m"

# Aggregate post-trim FastQC reports with MultiQC
if [ ! -s "${fastqc_trim_dir}/multiqc_report.html" ]; then
  "$conda_bin" run --no-capture-output -n multiqc "$multiqc" "$fastqc_trim_dir" --outdir "$fastqc_trim_dir"
fi
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc_trim: multiqc done\e[m"

# Seqkit Q30/Q20 stats on trimmed paired reads
if [ ! -s "${fastqc_trim_dir}/seqkit_stats.txt" ]; then
  "$conda_bin" run --no-capture-output -n seqkit "$seqkit" stats "${trim_dir}/"*_paired.fq.gz -aT -j "$MAX_Q30_STAT_THREADS" \
    -o "${fastqc_trim_dir}/seqkit_stats.txt"
fi
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc_trim: seqkit stats done\e[m"
