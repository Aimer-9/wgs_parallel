#!/bin/bash
# Stage 01: raw-read quality control.
#
# Called by wgs.sh; positional arguments:
#   1 work_dir                  Pipeline working directory.
#   2 sample_manifest_tsv       Normalized sample metadata from config.py.
#   3 fastqc                    FastQC executable or absolute path.
#   4 CPU fraction              Fraction of allocated CPUs available to stage.
#   5 max sample-pair jobs      Upper bound before R1/R2 expansion (default 10).
#   6 seqkit threads            Threads used for aggregate read statistics.
#   7 FastQC extra_args         Shell-style optional FastQC arguments.
#   8 multiqc command           Command executed in Conda env "wgs_parallel".
#   9 seqkit command            Command executed in Conda env "wgs_parallel".
#  10 conda executable          Conda used for MultiQC and SeqKit.
#
# Inputs are discovered from each sample's raw_dir in the manifest. R1 and R2
# are separate one-thread FastQC tasks because FastQC -t controls simultaneous
# files, not CPU threads per file. Outputs, the task list, MultiQC report, and
# SeqKit table are written to output/01_fastqc. Existing reports are resumable;
# checkpoint recovery removes only an interrupted file's HTML/ZIP pair.
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

# Materialize absolute FASTQ paths once so GNU Parallel workers do not repeat
# manifest lookup and so missing pairs fail before any FastQC job is launched.
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

# Aggregate every raw FastQC result into a single interactive HTML report.
if [ ! -s "${fastqc_dir}/multiqc_report.html" ]; then
  "$conda_bin" run --no-capture-output -n wgs_parallel "$multiqc" "$fastqc_dir" --outdir "$fastqc_dir"
fi
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc: multiqc done\e[m"

# Calculate one tabular SeqKit summary across the exact R1/R2 inputs above.
if [ ! -s "${fastqc_dir}/seqkit_stats.txt" ]; then
  seqkit_inputs=()
  while IFS= read -r sample_id; do
    [ -z "$sample_id" ] && continue
    seqkit_inputs+=("$(find_sample_read "$sample_id" "$sample_manifest_tsv" R1)")
    seqkit_inputs+=("$(find_sample_read "$sample_id" "$sample_manifest_tsv" R2)")
  done < "${work_dir}/sample_list.txt"
  "$conda_bin" run --no-capture-output -n wgs_parallel "$seqkit" stats "${seqkit_inputs[@]}" -aT -j "$MAX_Q30_STAT_THREADS" \
    -o "${fastqc_dir}/seqkit_stats.txt"
fi
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mfastqc: seqkit stats done\e[m"
