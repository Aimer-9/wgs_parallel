#!/bin/bash
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "${REPO_DIR}/scripts/00_util.sh"

work_dir=$1
hsa_version=$2
ref_fa=$3
bwamem=$4
samtools=$5
mosdepth_path=$6
MAX_PROCESSOR_USE_PERCENT=$7
MAX_BWA_JOBS=${8:-3}      # BWA-MEM2 is RAM-heavy; keep low
MAX_SORT_JOBS=${9:-10}
MAX_INDEX_JOBS=${10:-10}
MAX_MOSDEPTH_JOBS=${11:-10}
BWA_EXTRA_ARGS=${12:--v 1 -M}
SORT_MEMORY_PER_THREAD=${13:-2G}
SORT_EXTRA_ARGS=${14:-}
mosdepth=${15:-mosdepth}
MOSDEPTH_THREADS=${16:-4}
MOSDEPTH_EXTRA_ARGS=${17:--n}
conda_bin=${18:-conda}
SORT_THREADS_PER_JOB=${19:-2}

num_samples=$(wc -l < "${work_dir}/sample_list.txt")
nproc=$(nproc)
for var in MAX_BWA_JOBS MAX_SORT_JOBS MAX_INDEX_JOBS MAX_MOSDEPTH_JOBS; do
  [ "${!var}" -gt "$num_samples" ] && eval "$var=$num_samples"
done

threads_bwa=$(calculate_threads "$MAX_PROCESSOR_USE_PERCENT" "$MAX_BWA_JOBS")
threads_index=$(calculate_threads "$MAX_PROCESSOR_USE_PERCENT" "$MAX_INDEX_JOBS")
max_threads_sort=$(calculate_threads "$MAX_PROCESSOR_USE_PERCENT" "$MAX_SORT_JOBS")
threads_sort=$SORT_THREADS_PER_JOB
[ "$threads_sort" -gt "$max_threads_sort" ] && threads_sort=$max_threads_sort

# Keep concurrent sort jobs within the CPU budget after applying the fixed
# per-job thread cap.
sort_jobs_by_cpu=$(( $(calculate_cpu_budget "$MAX_PROCESSOR_USE_PERCENT") / threads_sort ))
[ "$sort_jobs_by_cpu" -lt 1 ] && sort_jobs_by_cpu=1
[ "$MAX_SORT_JOBS" -gt "$sort_jobs_by_cpu" ] && MAX_SORT_JOBS=$sort_jobs_by_cpu

bwa_dir="${work_dir}/output/04_bwa"
if [ ! -d "$bwa_dir" ]; then mkdir -p "$bwa_dir"; fi

# Build BWA-MEM2 index if missing (not included in GATK bundle)
if [ ! -s "${ref_fa}.bwt.2bit.64" ]; then
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mbwa index: building\e[m"
  "$bwamem" index "$ref_fa"
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mbwa index: done\e[m"
fi

# ── Align ────────────────────────────────────────────────────────────────────
run_bwa() {
  local sample_id=$1
  local trim_1="${work_dir}/output/02_trim/${sample_id}_1_paired.fq.gz"
  local trim_2="${work_dir}/output/02_trim/${sample_id}_2_paired.fq.gz"
  local out_bam="${bwa_dir}/${sample_id}.bam"
  local -a extra_args=()
  read_extra_args "$BWA_EXTRA_ARGS" extra_args || return 1

  if [ ! -s "$trim_1" ] || [ ! -s "$trim_2" ]; then
    echo "Error: trimmed FASTQ not found for ${sample_id} — run 02_trim.sh first"
    return 1
  fi
  if [ -s "$out_bam" ]; then
    echo -e "\e[37;42mbwa: ${sample_id} exists\e[m"; return 0
  fi

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mbwa: ${sample_id} processing\e[m"
  "$bwamem" mem "${extra_args[@]}" \
    -t "$threads_bwa" \
    -R "@RG\tID:lane\tPL:ILLUMINA\tLB:library\tSM:${sample_id}" \
    "$ref_fa" "$trim_1" "$trim_2" \
    | "$samtools" view -S -b > "$out_bam"
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mbwa: ${sample_id} done\e[m"
}

cleanup_run_bwa() {
  rm -f "${bwa_dir}/$1.bam"
}

date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mbwa: start\e[m"
parallel_run_sample_list "${work_dir}/sample_list.txt" "$MAX_BWA_JOBS" run_bwa bwa || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mbwa: all done\e[m"

# ── Sort ─────────────────────────────────────────────────────────────────────
run_sort() {
  local sample_id=$1
  local in_bam="${bwa_dir}/${sample_id}.bam"
  local out_bam="${bwa_dir}/${sample_id}.sorted.bam"
  local -a extra_args=()
  read_extra_args "$SORT_EXTRA_ARGS" extra_args || return 1

  if [ ! -s "$in_bam" ]; then
    echo "Error: ${in_bam} not found — bwa step may have failed"; return 1
  fi
  if [ -s "$out_bam" ]; then
    echo -e "\e[37;42msort: ${sample_id} exists\e[m"; return 0
  fi

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42msort: ${sample_id} processing\e[m"
  "$samtools" sort -@ "$threads_sort" -m "$SORT_MEMORY_PER_THREAD" -O bam "${extra_args[@]}" \
    -o "$out_bam" "$in_bam"
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42msort: ${sample_id} done\e[m"
}

cleanup_run_sort() {
  local out_bam="${bwa_dir}/$1.sorted.bam"
  rm -f "$out_bam" "${out_bam}.tmp."*
}

date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42msort: start\e[m"
echo "Sort resources: ${MAX_SORT_JOBS} parallel jobs, ${threads_sort} threads/job, ${SORT_MEMORY_PER_THREAD}/thread"
parallel_run_sample_list "${work_dir}/sample_list.txt" "$MAX_SORT_JOBS" run_sort sort || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42msort: all done\e[m"

# ── Index ────────────────────────────────────────────────────────────────────
run_index() {
  local sample_id=$1
  local sorted_bam="${bwa_dir}/${sample_id}.sorted.bam"

  if [ ! -s "$sorted_bam" ]; then
    echo "Error: ${sorted_bam} not found — sort step may have failed"; return 1
  fi
  if [ -s "${sorted_bam}.bai" ]; then
    echo -e "\e[37;42mindex: ${sample_id} exists\e[m"; return 0
  fi

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mindex: ${sample_id} processing\e[m"
  "$samtools" index -@ "$threads_index" "$sorted_bam"
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mindex: ${sample_id} done\e[m"
}

cleanup_run_index() {
  rm -f "${bwa_dir}/$1.sorted.bam.bai" "${bwa_dir}/$1.sorted.bai"
}

date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mindex: start\e[m"
parallel_run_sample_list "${work_dir}/sample_list.txt" "$MAX_INDEX_JOBS" run_index index || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mindex: all done\e[m"

# ── Mosdepth coverage ────────────────────────────────────────────────────────
run_mosdepth() {
  local sample_id=$1
  local sorted_bam="${bwa_dir}/${sample_id}.sorted.bam"
  local out_prefix="${bwa_dir}/${sample_id}.sorted"
  local -a extra_args=()
  read_extra_args "$MOSDEPTH_EXTRA_ARGS" extra_args || return 1

  if [ ! -s "$sorted_bam" ]; then
    echo "Error: ${sorted_bam} not found — sort step may have failed"; return 1
  fi
  if [ -s "${out_prefix}.mosdepth.global.dist.txt" ]; then
    echo -e "\e[37;42mmosdepth: ${sample_id} exists\e[m"; return 0
  fi

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mmosdepth: ${sample_id} processing\e[m"
  # mosdepth gains little above 4 threads
  "$conda_bin" run --no-capture-output -n mosdepth "$mosdepth" -t "$MOSDEPTH_THREADS" "${extra_args[@]}" "$out_prefix" "$sorted_bam"
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mmosdepth: ${sample_id} done\e[m"
}

cleanup_run_mosdepth() {
  local out_prefix="${bwa_dir}/$1.sorted"
  rm -f "${out_prefix}.mosdepth."* "${out_prefix}.per-base.bed.gz"* \
    "${out_prefix}.regions.bed.gz"* "${out_prefix}.thresholds.bed.gz"* \
    "${out_prefix}.quantized.bed.gz"*
}

date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mmosdepth: start\e[m"
parallel_run_sample_list "${work_dir}/sample_list.txt" "$MAX_MOSDEPTH_JOBS" run_mosdepth mosdepth || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mmosdepth: all done\e[m"

# Coverage comparison HTML across all samples
if [ ! -s "${bwa_dir}/compare.coverage_analysis.html" ]; then
  python3 "${mosdepth_path}/scripts/plot-dist.py" \
    -o "${bwa_dir}/compare.coverage_analysis.html" \
    "${bwa_dir}"/*.sorted.mosdepth.global.dist.txt \
    > "${bwa_dir}/compare.coverage_analysis.txt"
fi

# ── Coverage summary table (10x / 20x / 30x per sample) ─────────────────────
output_stat="${bwa_dir}/coverage_stat.txt"
if [ -s "$output_stat" ]; then
  echo -e "\e[37;42mcoverage stat: exists\e[m"
else
  echo -e "Sample_ID\tMean_Depth\tCoverage_10x(%)\tCoverage_20x(%)\tCoverage_30x(%)" > "$output_stat"
  while IFS= read -r sample_id; do
    local_dist="${bwa_dir}/${sample_id}.sorted.mosdepth.global.dist.txt"
    local_summary="${bwa_dir}/${sample_id}.sorted.mosdepth.summary.txt"
    if [ ! -s "$local_dist" ]; then
      echo "Error: mosdepth output not found for ${sample_id}"; continue
    fi
    mean_depth=$(awk -F'\t' '$1=="total"{print $4}' "$local_summary")
    cov_10x=$(awk -F'\t' '$1=="total" && $2==10{printf "%.2f", 100*$3}' "$local_dist")
    cov_20x=$(awk -F'\t' '$1=="total" && $2==20{printf "%.2f", 100*$3}' "$local_dist")
    cov_30x=$(awk -F'\t' '$1=="total" && $2==30{printf "%.2f", 100*$3}' "$local_dist")
    echo -e "${sample_id}\t${mean_depth}\t${cov_10x}\t${cov_20x}\t${cov_30x}" >> "$output_stat"
  done < "${work_dir}/sample_list.txt"
fi
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mcoverage stat: done\e[m"
