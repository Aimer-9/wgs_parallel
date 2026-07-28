#!/bin/bash

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPTS_DIR="${REPO_DIR}/scripts"

usage() {
  cat <<'EOF'
Usage: bash wgs.sh <command> [sample.csv] [config.yaml]

Defaults:
  sample.csv   -> example/config/sample.csv
  config.yaml  -> ./config.yaml

Commands:
  qc     QC + trimming, then stop for report review
  run    Alignment -> variant calling -> annotation (requires QC outputs)
  all    Run QC, then continue into alignment -> annotation
  stop   Stop all running pipeline jobs started by this workspace
EOF
}

set_input_paths() {
  if [ $# -gt 2 ]; then
    usage
    exit 1
  fi

  if [ $# -ge 1 ]; then
    export SAMPLE_CSV=$1
  fi
  if [ $# -ge 2 ]; then
    export CONFIG_YAML=$2
  fi
}

conda_env_exists() {
  local env_name=$1
  "${conda_bin:-conda}" env list 2>/dev/null | awk '{print $1}' | grep -qx "$env_name"
}

conda_env_executable_exists() {
  local env_name=$1
  local executable=$2
  "${conda_bin:-conda}" run --no-capture-output -n "$env_name" \
    bash -c 'command -v "$1" >/dev/null 2>&1' _ "$executable" >/dev/null 2>&1
}

executable_exists() {
  local value=$1
  if [[ "$value" == */* ]]; then
    [ -x "$value" ]
  else
    command -v "$value" >/dev/null 2>&1
  fi
}

validate_qc_config() {
  local errors=()

  for var in work_dir sample_manifest_tsv fastqc trim_path MAX_PROCESSOR_USE_PERCENT; do
    if [ -z "${!var:-}" ]; then
      errors+=("  [params] '$var' is not set")
    fi
  done

  if [ -n "${sample_manifest_tsv:-}" ] && [ ! -f "$sample_manifest_tsv" ]; then
    errors+=("  [paths] sample manifest not found: $sample_manifest_tsv")
  fi
  if [ -n "${fastqc:-}" ] && ! executable_exists "$fastqc"; then
    errors+=("  [tools] fastqc not found or not executable: $fastqc")
  fi
  if [ -n "${trim_path:-}" ] && [ ! -f "${trim_path}/trimmomatic-0.39.jar" ]; then
    errors+=("  [tools] Trimmomatic jar not found: ${trim_path}/trimmomatic-0.39.jar")
  fi

  for tool_var in java conda_bin parallel_bin; do
    local tool_path="${!tool_var:-}"
    if [ -z "$tool_path" ] || ! executable_exists "$tool_path"; then
      errors+=("  [tools] '$tool_var' not found or not executable: $tool_path")
    fi
  done
  if [ -n "${conda_bin:-}" ] && executable_exists "$conda_bin"; then
    if ! conda_env_executable_exists multiqc "${multiqc:-multiqc}"; then
      errors+=("  [conda] '${multiqc:-multiqc}' not found in environment: multiqc")
    fi
    if ! conda_env_executable_exists seqkit "${seqkit:-seqkit}"; then
      errors+=("  [conda] '${seqkit:-seqkit}' not found in environment: seqkit")
    fi
  fi

  if [ ${#errors[@]} -gt 0 ]; then
    echo -e "\e[37;41mConfiguration errors - fix config.yaml before re-running:\e[m"
    for msg in "${errors[@]}"; do
      echo "$msg"
    done
    exit 1
  fi
}

validate_run_config() {
  local errors=()

  for var in work_dir sample_manifest_tsv \
             bwamem samtools gatk bcftools annovar_path mosdepth_path \
             hsa_version MAX_PROCESSOR_USE_PERCENT update_clinvar \
             ref_fa known_indels dbsnp Mills hapmap file_1000G_omni file_1000G_phase1; do
    if [ -z "${!var:-}" ]; then
      errors+=("  [params] '$var' is not set")
    fi
  done

  if [ -n "${sample_manifest_tsv:-}" ] && [ ! -f "$sample_manifest_tsv" ]; then
    errors+=("  [paths] sample manifest not found: $sample_manifest_tsv")
  fi
  if [ ! -d "$SCRIPTS_DIR" ]; then
    errors+=("  [paths] scripts/ directory not found: $SCRIPTS_DIR")
  fi
  if [ -n "${annovar_path:-}" ] && [ ! -d "$annovar_path" ]; then
    errors+=("  [paths] annovar_path directory not found: $annovar_path")
  fi
  if [ -n "${work_dir:-}" ]; then
    local work_parent
    work_parent="$(dirname "$work_dir")"
    if [ ! -d "$work_parent" ]; then
      errors+=("  [paths] parent of work_dir does not exist: $work_parent")
    elif [ ! -w "$work_parent" ]; then
      errors+=("  [paths] parent of work_dir is not writable: $work_parent")
    fi
  fi

  for tool_var in bwamem samtools gatk bcftools; do
    local tool_path="${!tool_var:-}"
    if [ -n "$tool_path" ] && ! executable_exists "$tool_path"; then
      errors+=("  [tools] '$tool_var' not found or not executable: $tool_path")
    fi
  done
  if [ -n "${mosdepth_path:-}" ] && [ ! -f "${mosdepth_path}/scripts/plot-dist.py" ]; then
    errors+=("  [tools] mosdepth plot-dist.py not found: ${mosdepth_path}/scripts/plot-dist.py")
  fi

  for tool_var in java conda_bin parallel_bin; do
    local tool_path="${!tool_var:-}"
    if [ -z "$tool_path" ] || ! executable_exists "$tool_path"; then
      errors+=("  [tools] '$tool_var' not found or not executable: $tool_path")
    fi
  done
  if [ -n "${conda_bin:-}" ] && executable_exists "$conda_bin"; then
    if ! conda_env_executable_exists mosdepth "${mosdepth:-mosdepth}"; then
      errors+=("  [conda] '${mosdepth:-mosdepth}' not found in environment: mosdepth")
    fi
  fi
  if [ "${update_clinvar:-no}" = "yes" ] && ! command -v conda >/dev/null 2>&1; then
    errors+=("  [conda] conda is required when update_clinvar=yes")
  fi

  for ref_var in ref_fa known_indels dbsnp Mills hapmap file_1000G_omni file_1000G_phase1; do
    local ref_path="${!ref_var:-}"
    if [ -n "$ref_path" ] && [ ! -f "$ref_path" ]; then
      errors+=("  [refs] '$ref_var' not found: $ref_path")
    fi
  done

  if awk -F '\t' 'NF >= 6 && $6 != "" { found=1 } END { exit !found }' "$sample_manifest_tsv"; then
    for ref_var in germline_resource panel_of_normals; do
      local ref_path="${!ref_var:-}"
      if [ -z "$ref_path" ]; then
        errors+=("  [somatic] '$ref_var' is not set (required by normal_bam rows)")
      elif [ ! -f "$ref_path" ]; then
        errors+=("  [somatic] '$ref_var' not found: $ref_path")
      fi
    done
  fi

  if [ ${#errors[@]} -gt 0 ]; then
    echo -e "\e[37;41mConfiguration errors - fix config.yaml before re-running:\e[m"
    for msg in "${errors[@]}"; do
      echo "$msg"
    done
    exit 1
  fi

  echo -e "\e[37;42mValidation passed\e[m"
}

prepare_workspace() {
  prepare_sample_inputs || exit 1
  echo "Total number of samples: $(wc -l < "${work_dir}/sample_list.txt")"
}

ensure_trimmed_fastqs() {
  local missing=()

  while IFS= read -r sample_id; do
    [ -z "$sample_id" ] && continue
    if [ ! -s "${work_dir}/output/02_trim/${sample_id}_1_paired.fq.gz" ] || [ ! -s "${work_dir}/output/02_trim/${sample_id}_2_paired.fq.gz" ]; then
      missing+=("$sample_id")
    fi
  done < "${work_dir}/sample_list.txt"

  if [ ${#missing[@]} -gt 0 ]; then
    echo -e "\e[37;41mError: trimmed FASTQs not found in ${work_dir}/output/02_trim for:\e[m"
    printf '  %s\n' "${missing[@]}"
    echo "Run: bash wgs.sh qc ${SAMPLE_CSV:-example/config/sample.csv} ${CONFIG_YAML:-config.yaml}"
    exit 1
  fi
}

run_qc_steps() {
  bash "${SCRIPTS_DIR}/01_fastqc.sh"      "$work_dir" "$sample_manifest_tsv" "$fastqc" "$MAX_PROCESSOR_USE_PERCENT" "$MAX_FASTQC_RAW_JOBS" "$FASTQC_RAW_Q30_THREADS" "$FASTQC_RAW_EXTRA_ARGS" "$multiqc" "$seqkit" "$conda_bin" || return 1
  bash "${SCRIPTS_DIR}/02_trim.sh"        "$work_dir" "$sample_manifest_tsv" "$trim_path" "$MAX_PROCESSOR_USE_PERCENT" "$MAX_TRIM_JOBS" "$TRIM_LEADING" "$TRIM_TRAILING" "$TRIM_MINLEN" "$TRIM_SLIDINGWINDOW" "$TRIM_EXTRA_ARGS" "$java" || return 1
  bash "${SCRIPTS_DIR}/03_fastqc_trim.sh" "$work_dir" "$fastqc" "$MAX_PROCESSOR_USE_PERCENT" "$MAX_FASTQC_TRIM_JOBS" "$FASTQC_TRIM_Q30_THREADS" "$FASTQC_TRIM_EXTRA_ARGS" "$multiqc" "$seqkit" "$conda_bin" || return 1

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mQC complete\e[m"
  echo ""
  echo "Review the following reports before variant calling:"
  echo "  Raw FastQC / MultiQC  : ${work_dir}/output/01_fastqc/multiqc_report.html"
  echo "  Trimmed FastQC / MultiQC: ${work_dir}/output/03_fastqc_trim/multiqc_report.html"
  echo "  Raw seqkit stats      : ${work_dir}/output/01_fastqc/seqkit_stats.txt"
  echo "  Trimmed seqkit stats  : ${work_dir}/output/03_fastqc_trim/seqkit_stats.txt"
}

run_variant_steps() {
  mkdir -p "${work_dir}/output"
  prepare_primary_contigs || return 1

  bash "${SCRIPTS_DIR}/04_bwa.sh" "$work_dir" "$hsa_version" "$ref_fa" "$bwamem" "$samtools" "$mosdepth_path" "$MAX_PROCESSOR_USE_PERCENT" \
    "$MAX_BWA_JOBS" "$MAX_SORT_JOBS" "$MAX_INDEX_JOBS" "$MAX_MOSDEPTH_JOBS" "$BWA_EXTRA_ARGS" "$SORT_MEMORY_PER_THREAD" "$SORT_EXTRA_ARGS" "$mosdepth" "$MOSDEPTH_THREADS" "$MOSDEPTH_EXTRA_ARGS" "$conda_bin" "$SORT_THREADS_PER_JOB" || return 1
  bash "${SCRIPTS_DIR}/05_markduplicated.sh" "$work_dir" "$hsa_version" "$gatk" "$samtools" "$MAX_PROCESSOR_USE_PERCENT" "$MAX_MARKDUP_JOBS" "$MARKDUP_EXTRA_ARGS" || return 1
  bash "${SCRIPTS_DIR}/06_bqsr.sh" "$work_dir" "$hsa_version" "$gatk" "$samtools" "$ref_fa" "$known_indels" "$dbsnp" "$Mills" "$MAX_PROCESSOR_USE_PERCENT" \
    "$MAX_BQSR_JOBS" "$BQSR_MIN_OUTPUT_BYTES" "$BQSR_EXTRA_ARGS" || return 1

  bash "${SCRIPTS_DIR}/07_HaplotypeCaller.sh" "$work_dir" "$hsa_version" "$gatk" "$ref_fa" "$MAX_HC_JOBS" "$primary_contigs_tsv" "$HC_NATIVE_THREADS" "$HC_KEEP_SCATTER" "$HC_EXTRA_ARGS" || return 1
  bash "${SCRIPTS_DIR}/08_vqsr.sh"            "$work_dir" "$hsa_version" "$gatk" "$ref_fa" \
    "$known_indels" "$dbsnp" "$Mills" "$hapmap" "$file_1000G_omni" "$file_1000G_phase1" \
    "$MAX_VQSR_JOBS" "$MAX_GENOTYPE_JOBS" "$GENOTYPE_EXTRA_ARGS" "$VQSR_SNP_TRUTH_SENSITIVITY" "$VQSR_INDEL_TRUTH_SENSITIVITY" "$VQSR_EXTRA_ARGS" "$bcftools" || return 1

  bash "${SCRIPTS_DIR}/09_annovar.sh" \
    "$work_dir" "$hsa_version" "$annovar_path" \
    "${work_dir}/output/08_VQSR" \
    ".indel.SNP.recalibrated.PASS.vcf" \
    "$MAX_PROCESSOR_USE_PERCENT" "$update_clinvar" "$MAX_ANNOVAR_JOBS" "${work_dir}/sample_list.txt" "$ANNOVAR_PROTOCOL" "$ANNOVAR_OPERATION" "$ANNOVAR_EXTRA_ARGS" \
    "${work_dir}/output/09_annovar" || return 1

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mWGS germline pipeline: all done\e[m"

  if awk -F '\t' 'NF >= 6 && $6 != "" { found=1 } END { exit !found }' "$sample_manifest_tsv"; then
    local mutect_status=0
    echo "normal_bam entries found - proceeding with Mutect2 tumor-normal calling"
    bash "${SCRIPTS_DIR}/10_Mutect2.sh" \
      "$work_dir" "$sample_manifest_tsv" "$gatk" "$samtools" "$bcftools" "$hsa_version" "$ref_fa" \
      "$germline_resource" "$panel_of_normals" "$primary_contigs_tsv" "$MAX_MUTECT2_JOBS" "$MUTECT2_EXTRA_ARGS" "$FILTER_MUTECT_EXTRA_ARGS" || mutect_status=$?
    bash "${SCRIPTS_DIR}/09_annovar.sh" \
      "$work_dir" "$hsa_version" "$annovar_path" \
      "${work_dir}/output/10_Mutect2" \
      ".mutect2.filtered.PASS.vcf" \
      "$MAX_PROCESSOR_USE_PERCENT" "$update_clinvar" "$MAX_ANNOVAR_JOBS" "${work_dir}/output/10_Mutect2/success_samples.txt" "$ANNOVAR_PROTOCOL" "$ANNOVAR_OPERATION" "$ANNOVAR_EXTRA_ARGS" \
      "${work_dir}/output/11_annovar_somatic" || return 1
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mWGS somatic pipeline: all done\e[m"
    [ "$mutect_status" -eq 0 ] || return "$mutect_status"
  else
    echo "No normal_bam entries in sample.csv - skipping Mutect2"
  fi
}

main() {
  local command=${1:-}
  if [ -z "$command" ]; then
    usage
    exit 1
  fi
  shift

  case "$command" in
    stop)
      source "${REPO_DIR}/scripts/00_util.sh"
      stop_pipeline_process_group "${REPO_DIR}"
      cleanup_pipeline_control_files "${REPO_DIR}"
      exit 0
      ;;
    qc|run|all) ;;
    *)
      usage
      exit 1
      ;;
  esac

  set_input_paths "$@"

  source "${REPO_DIR}/scripts/00_util.sh"
  load_pipeline_config "${REPO_DIR}" || exit 1
  init_run_logging "$command" || exit 1
  init_pipeline_control_files "${REPO_DIR}" || exit 1
  trap 'echo "Interrupt received, terminating pipeline..."; stop_pipeline_process_group "${REPO_DIR}"; cleanup_pipeline_control_files "${REPO_DIR}"; exit 130' INT TERM
  trap 'cleanup_pipeline_control_files "${REPO_DIR}"' EXIT

  case "$command" in
    qc)
      validate_qc_config
      prepare_workspace
      run_qc_steps
      ;;
    run)
      validate_run_config
      prepare_workspace
      ensure_trimmed_fastqs
      run_variant_steps
      ;;
    all)
      validate_qc_config
      validate_run_config
      prepare_workspace
      run_qc_steps
      ensure_trimmed_fastqs
      run_variant_steps
      ;;
  esac
}

main "$@"
