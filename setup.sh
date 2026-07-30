#!/bin/bash
# Bootstrap and validate the software/reference environment used by wgs.sh.
#
# The default setup performs three independently selectable operations:
#   1. install command-line tools beneath ./tools;
#   2. create the reproducible wgs_parallel Conda environment containing Python,
#      MultiQC, SeqKit, mosdepth, and ClinVar updater dependencies;
#   3. download hg19 and/or hg38 references beneath ./ref.
#
# Machine-detected paths are merged with example/config/config.yaml to produce
# the user's config.yaml without discarding documented parameter defaults.
# Downloads and installations are restart-aware and validation can be run alone
# with no changes. Interactive permission prompts precede optional downloads;
# software tasks are queued and may execute concurrently where safe.
#
# Common examples:
#   bash setup.sh
#   bash setup.sh download-tools --tools-dir /opt/wgs-tools
#   bash setup.sh download-ref --build 38 --ref-dir /data/references
#   bash setup.sh validate-only --validate-tools

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
PREFIX="${REPO_DIR}/tools"
REFS_DIR="${REPO_DIR}/ref"
BUILD="all"

DO_SOFTWARE=true
DO_CONDA=true
DO_REFS=true
VALIDATE_ONLY=false
VALIDATE_TOOLS=true
VALIDATE_REFS=true
SUBCOMMAND="setup"

usage() {
  cat <<'EOF'
Usage: bash setup.sh [subcommand] [options]

Subcommands:
  validate-only      validate existing tools and/or references, no install/download
  download-tools     download/install software tools only
  download-ref       download reference files only
  (none)             run full setup (software + conda env + refs)

Options:
  --tools-dir <dir>  tools install root (alias: --prefix, default: ./tools)
  --ref-dir <dir>    reference data root (alias: --refs-dir, default: ./ref)
  --build <19|38|all> reference builds to configure/download/validate (default: all)

Legacy mode switches:
  --software-only      install binaries only
  --conda-only         create conda environments only
  --refs-only          download reference files only

Validate-only selectors:
  --validate-tools   validate tools only
  --validate-refs    validate references only
EOF
}

if [ $# -gt 0 ]; then
  case "$1" in
    validate-only|download-tools|download-ref)
      SUBCOMMAND=$1
      shift
      ;;
  esac
fi

case "$SUBCOMMAND" in
  validate-only)
    DO_SOFTWARE=false
    DO_CONDA=false
    DO_REFS=false
    VALIDATE_ONLY=true
    ;;
  download-tools)
    DO_SOFTWARE=true
    DO_CONDA=false
    DO_REFS=false
    ;;
  download-ref)
    DO_SOFTWARE=false
    DO_CONDA=false
    DO_REFS=true
    ;;
  setup) ;;
  *)
    echo "Error: unknown subcommand: $SUBCOMMAND"
    usage
    exit 1
    ;;
esac

while [ $# -gt 0 ]; do
  case "$1" in
    --software-only)
      DO_CONDA=false
      DO_REFS=false
      ;;
    --conda-only)
      DO_SOFTWARE=false
      DO_REFS=false
      ;;
    --refs-only)
      DO_SOFTWARE=false
      DO_CONDA=false
      ;;
    --prefix|--tools-dir)
      shift
      [ $# -gt 0 ] || { echo "Error: --prefix/--tools-dir requires a directory"; exit 1; }
      PREFIX=$1
      ;;
    --refs-dir|--ref-dir)
      shift
      [ $# -gt 0 ] || { echo "Error: --refs-dir/--ref-dir requires a directory"; exit 1; }
      REFS_DIR=$1
      ;;
    --build)
      shift
      [ $# -gt 0 ] || { echo "Error: --build requires 19, 38, or all"; exit 1; }
      BUILD=$1
      ;;
    --validate-tools)
      VALIDATE_TOOLS=true
      VALIDATE_REFS=false
      ;;
    --validate-refs)
      VALIDATE_TOOLS=false
      VALIDATE_REFS=true
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      usage
      exit 1
      ;;
  esac
  shift
done

if ! $VALIDATE_ONLY && { ! $VALIDATE_TOOLS || ! $VALIDATE_REFS; }; then
  echo "Error: --validate-tools/--validate-refs are only supported with 'validate-only'."
  exit 1
fi

case "$BUILD" in
  19|38|all) ;;
  *)
    echo "Error: --build must be 19, 38, or all (got: $BUILD)"
    exit 1
    ;;
esac

declare -a SELECTED_BUILDS=()
case "$BUILD" in
  19) SELECTED_BUILDS=("19") ;;
  38) SELECTED_BUILDS=("38") ;;
  all) SELECTED_BUILDS=("19" "38") ;;
esac

fastqc="${PREFIX}/FastQC/fastqc"
trim_path="${PREFIX}/Trimmomatic-0.39"
bwamem="${PREFIX}/bwa-mem2-2.3_x64-linux/bwa-mem2"
samtools="${PREFIX}/samtools-1.22.1/samtools"
gatk="${PREFIX}/gatk-4.6.2.0/gatk"
annovar_path="${PREFIX}/annovar"
mosdepth_path="${PREFIX}/mosdepth"

refs_dir_19="${REFS_DIR}/gcp-public-data--broad-references-hg19/v0"
somatic_dir_19="${refs_dir_19}/somatic_hg19"
refs_dir_38="${REFS_DIR}/gcp-public-data--broad-references-hg38/v0"
somatic_dir_38="${refs_dir_38}/somatic_hg38"

PASS="\e[37;42m OK \e[m"
SKIP="\e[37;44m SKIP \e[m"
FAIL="\e[37;41m FAIL \e[m"

log() {
  date +"%Y-%m-%d %H:%M:%S"
  echo -e "$*"
}

need() {
  if ! command -v "$1" &>/dev/null; then
    echo -e "${FAIL}  required command not found: $1${2:+ - install $2}"
    exit 1
  fi
}

need_parallel() {
  need parallel "GNU parallel"
}

is_yes() {
  local ans=${1:-}
  [[ "$ans" =~ ^[Yy]([Ee][Ss])?$ ]]
}

ask_permission() {
  local label=$1
  local answer=""
  if [ ! -t 0 ]; then
    echo -e "${SKIP}  non-interactive shell: skip prompt for ${label}"
    return 1
  fi
  read -r -p "Install ${label}? [y/N] " answer
  is_yes "$answer"
}

download() {
  local url=$1
  local dest=$2
  if [ -f "$dest" ]; then
    return 0
  fi
  mkdir -p "$(dirname "$dest")"
  wget -q --show-progress -O "$dest" "$url" || {
    rm -f "$dest"
    return 1
  }
}

conda_env_exists() {
  local env_name=$1
  conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$env_name"
}

ensure_install_roots() {
  mkdir -p "$PREFIX" "$REFS_DIR" "$annovar_path" "$mosdepth_path"
}

parallel_eval_tab_jobs() {
  local jobs_file=$1
  local max_jobs=${2:-0}
  local v=""
  local f=""
  if [ ! -s "$jobs_file" ]; then
    return 0
  fi
  local jobs_opt=()
  if [ "$max_jobs" -gt 0 ] 2>/dev/null; then
    jobs_opt=(--jobs "$max_jobs")
  fi

  while IFS= read -r v; do
    case "$v" in
      BASHOPTS|BASHPID|BASH_ALIASES|BASH_ARGC|BASH_ARGV|BASH_CMDS|BASH_COMMAND|BASH_LINENO|BASH_SOURCE|BASH_VERSINFO|DIRSTACK|EUID|GROUPS|FUNCNAME|LINENO|PIPESTATUS|PPID|SHELLOPTS|UID)
        continue
        ;;
    esac
    export "$v" 2>/dev/null || true
  done < <(compgen -v)

  while IFS= read -r f; do
    export -f "$f" 2>/dev/null || true
  done < <(declare -F | awk '{print $3}')

  parallel --will-cite \
    --colsep '\t' \
    --line-buffer \
    --halt soon,fail=1 \
    "${jobs_opt[@]}" \
    bash -lc 'echo "Running: $1"; eval "$2" && echo -e "${PASS}  $1" || { echo -e "${FAIL}  $1"; exit 1; }' _ {1} {2} \
    :::: "$jobs_file"
}

yaml_value_if_file() {
  local path=$1
  if [ -s "$path" ]; then
    echo "$path"
  else
    echo "null"
  fi
}

yaml_value_if_exec() {
  local path=$1
  if [ -x "$path" ]; then
    echo "$path"
  else
    echo "null"
  fi
}

yaml_value_if_exists() {
  local path=$1
  if [ -e "$path" ]; then
    echo "$path"
  else
    echo "null"
  fi
}

write_build_yaml_validated() {
  local build=$1

  if [ "$build" = "19" ]; then
    cat <<EOF
  "19":
    ref_fa: $(yaml_value_if_file "${refs_dir_19}/Homo_sapiens_assembly19.fasta")
    known_indels: $(yaml_value_if_file "${refs_dir_19}/Homo_sapiens_assembly19.known_indels_20120518.vcf")
    dbsnp: $(yaml_value_if_file "${refs_dir_19}/dbsnp_138.b37.vcf.gz")
    Mills: $(yaml_value_if_file "${refs_dir_19}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz")
    hapmap: $(yaml_value_if_file "${refs_dir_19}/hapmap_3.3.b37.vcf.gz")
    file_1000G_omni: $(yaml_value_if_file "${refs_dir_19}/1000G_omni2.5.b37.vcf.gz")
    file_1000G_phase1: $(yaml_value_if_file "${refs_dir_19}/1000G_phase1.snps.high_confidence.b37.vcf.gz")
    germline_resource: $(yaml_value_if_file "${somatic_dir_19}/af-only-gnomad.raw.sites.vcf.gz")
    panel_of_normals: $(yaml_value_if_file "${somatic_dir_19}/Mutect2-WGS-panel-b37.vcf.gz")
EOF
    return
  fi

  cat <<EOF
  "38":
    ref_fa: $(yaml_value_if_file "${refs_dir_38}/Homo_sapiens_assembly38.fasta")
    known_indels: $(yaml_value_if_file "${refs_dir_38}/Homo_sapiens_assembly38.known_indels.vcf.gz")
    dbsnp: $(yaml_value_if_file "${refs_dir_38}/Homo_sapiens_assembly38.dbsnp138.vcf.gz")
    Mills: $(yaml_value_if_file "${refs_dir_38}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
    hapmap: $(yaml_value_if_file "${refs_dir_38}/hapmap_3.3.hg38.vcf.gz")
    file_1000G_omni: $(yaml_value_if_file "${refs_dir_38}/1000G_omni2.5.hg38.vcf.gz")
    file_1000G_phase1: $(yaml_value_if_file "${refs_dir_38}/1000G_phase1.snps.high_confidence.hg38.vcf.gz")
    germline_resource: $(yaml_value_if_file "${somatic_dir_38}/af-only-gnomad.hg38.vcf.gz")
    panel_of_normals: $(yaml_value_if_file "${somatic_dir_38}/1000g_pon.hg38.vcf.gz")
EOF
}

write_config() {
  local out="${REPO_DIR}/config.yaml"
  local detected="${TMPDIR:-/tmp}/wgs_parallel.detected.$$.yaml"
  local fastqc_valid trim_valid bwa_valid samtools_valid gatk_valid annovar_valid mosdepth_valid
  local java_valid bcftools_valid parallel_valid multiqc_valid seqkit_valid mosdepth_exec

  fastqc_valid=$(yaml_value_if_exec "${fastqc}")
  trim_valid=$(yaml_value_if_file "${trim_path}/trimmomatic-0.39.jar")
  bwa_valid=$(yaml_value_if_exec "${bwamem}")
  samtools_valid=$(yaml_value_if_exec "${samtools}")
  gatk_valid=$(yaml_value_if_exec "${gatk}")
  annovar_valid=$(yaml_value_if_exists "${annovar_path}")
  mosdepth_valid=$(yaml_value_if_exists "${mosdepth_path}")
  java_valid=$(command -v java 2>/dev/null || echo java)
  local conda_valid
  conda_valid=$(command -v conda 2>/dev/null || echo conda)
  bcftools_valid=$(command -v bcftools 2>/dev/null || echo bcftools)
  parallel_valid=$(command -v parallel 2>/dev/null || echo parallel)
  multiqc_valid=$(command -v multiqc 2>/dev/null || true)
  seqkit_valid=$(command -v seqkit 2>/dev/null || true)
  mosdepth_exec=$(command -v mosdepth 2>/dev/null || true)
  if command -v conda >/dev/null 2>&1; then
    local conda_base
    conda_base=$(conda info --base 2>/dev/null || true)
    [ -z "$multiqc_valid" ] && multiqc_valid="${conda_base}/envs/wgs_parallel/bin/multiqc"
    [ -z "$seqkit_valid" ] && seqkit_valid="${conda_base}/envs/wgs_parallel/bin/seqkit"
    [ -z "$mosdepth_exec" ] && mosdepth_exec="${conda_base}/envs/wgs_parallel/bin/mosdepth"
  fi
  [ -n "$multiqc_valid" ] || multiqc_valid=multiqc
  [ -n "$seqkit_valid" ] || seqkit_valid=seqkit
  [ -n "$mosdepth_exec" ] || mosdepth_exec=mosdepth

  cat > "$detected" <<EOF
tools:
  java: ${java_valid}
  conda: ${conda_valid}
  fastqc: ${fastqc_valid}
  trimmomatic_dir: $(yaml_value_if_exists "${trim_path}")
  bwa_mem2: ${bwa_valid}
  samtools: ${samtools_valid}
  gatk: ${gatk_valid}
  bcftools: ${bcftools_valid}
  parallel: ${parallel_valid}
  multiqc: ${multiqc_valid}
  seqkit: ${seqkit_valid}
  mosdepth: ${mosdepth_exec}
  annovar_dir: ${annovar_valid}
  mosdepth_dir: ${mosdepth_valid}
references:
EOF

  for build in "${SELECTED_BUILDS[@]}"; do
    write_build_yaml_validated "$build" | sed \
      -e 's/^  "19":/  "19":/' -e 's/^  "38":/  "38":/' \
      -e 's/    ref_fa:/    fasta:/' \
      -e 's/    Mills:/    mills:/' \
      -e 's/    file_1000G_omni:/    omni_1000g:/' \
      -e 's/    file_1000G_phase1:/    phase1_1000g:/' >> "$detected"
  done

  python3 "${REPO_DIR}/scripts/config.py" merge-machine \
    --detected "$detected" \
    --template "${REPO_DIR}/example/config/config.yaml" \
    --output "$out"
  rm -f "$detected"
}

download_refs_for_build() {
  local build=$1
  local base=""
  local somatic=""
  local ref_dir=""
  local somatic_dir=""

  gsutil_get() {
    local gcs_path=$1
    local local_dir=$2
    local fname
    fname=$(basename "$gcs_path")
    if compgen -G "${local_dir}/${fname}*" > /dev/null; then
      echo -e "${SKIP}  $fname"
    else
      log "Downloading $fname..."
      gsutil -m cp "${gcs_path}*" "$local_dir/" \
        && echo -e "${PASS}  $fname" \
        || echo -e "${FAIL}  $fname"
    fi
  }

  if [ "$build" = "19" ]; then
    base="gs://gcp-public-data--broad-references/hg19/v0"
    somatic="gs://gatk-best-practices/somatic-b37"
    ref_dir="$refs_dir_19"
    somatic_dir="$somatic_dir_19"
  else
    base="gs://gcp-public-data--broad-references/hg38/v0"
    somatic="gs://gatk-best-practices/somatic-hg38"
    ref_dir="$refs_dir_38"
    somatic_dir="$somatic_dir_38"
  fi

  mkdir -p "$ref_dir" "$somatic_dir"
  log "Downloading reference files for hg${build}..."

  if [ "$build" = "19" ]; then
    gsutil_get "${base}/Homo_sapiens_assembly19.fasta"                     "$ref_dir"
    gsutil_get "${base}/Homo_sapiens_assembly19.known_indels_20120518.vcf" "$ref_dir"
    gsutil_get "${base}/dbsnp_138.b37.vcf.gz"                              "$ref_dir"
    gsutil_get "${base}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"   "$ref_dir"
    gsutil_get "${base}/hapmap_3.3.b37.vcf.gz"                             "$ref_dir"
    gsutil_get "${base}/1000G_omni2.5.b37.vcf.gz"                          "$ref_dir"
    gsutil_get "${base}/1000G_phase1.snps.high_confidence.b37.vcf.gz"      "$ref_dir"
    gsutil_get "${somatic}/af-only-gnomad.raw.sites.vcf.gz"                "$somatic_dir"
    gsutil_get "${somatic}/Mutect2-WGS-panel-b37.vcf.gz"                   "$somatic_dir"
  else
    gsutil_get "${base}/Homo_sapiens_assembly38.fasta"                     "$ref_dir"
    gsutil_get "${base}/Homo_sapiens_assembly38.known_indels.vcf.gz"       "$ref_dir"
    gsutil_get "${base}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"           "$ref_dir"
    gsutil_get "${base}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"  "$ref_dir"
    gsutil_get "${base}/hapmap_3.3.hg38.vcf.gz"                            "$ref_dir"
    gsutil_get "${base}/1000G_omni2.5.hg38.vcf.gz"                         "$ref_dir"
    gsutil_get "${base}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"     "$ref_dir"
    gsutil_get "${somatic}/af-only-gnomad.hg38.vcf.gz"                     "$somatic_dir"
    gsutil_get "${somatic}/1000g_pon.hg38.vcf.gz"                          "$somatic_dir"
  fi

  if [ -f "${annovar_path}/annotate_variation.pl" ]; then
    local humandb="${annovar_path}/humandb"
    mkdir -p "$humandb"
    log "Pre-downloading ANNOVAR databases for hg${build}..."
    for db in refGene cytoBand exac03 avsnp150; do
      if [ ! -s "${humandb}/hg${build}_${db}.txt" ]; then
        log "Downloading ANNOVAR db: $db"
        perl "${annovar_path}/annotate_variation.pl" \
          -buildver "hg${build}" -downdb -webfrom annovar "$db" "$humandb"
      else
        echo -e "${SKIP}  ANNOVAR db: $db (hg${build})"
      fi
    done
  else
    echo -e "\e[33m[SKIP]\e[m   ANNOVAR not installed - databases will be downloaded on first annotation run"
  fi
}

print_ref_summary_for_build() {
  local build=$1
  local ref_dir=""
  local somatic_dir=""
  local -a paths=()

  if [ "$build" = "19" ]; then
    ref_dir="$refs_dir_19"
    somatic_dir="$somatic_dir_19"
    paths=(
      "${ref_dir}/Homo_sapiens_assembly19.fasta"
      "${ref_dir}/Homo_sapiens_assembly19.known_indels_20120518.vcf"
      "${ref_dir}/dbsnp_138.b37.vcf.gz"
      "${ref_dir}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
      "${ref_dir}/hapmap_3.3.b37.vcf.gz"
      "${ref_dir}/1000G_omni2.5.b37.vcf.gz"
      "${ref_dir}/1000G_phase1.snps.high_confidence.b37.vcf.gz"
      "${somatic_dir}/af-only-gnomad.raw.sites.vcf.gz"
      "${somatic_dir}/Mutect2-WGS-panel-b37.vcf.gz"
    )
  else
    ref_dir="$refs_dir_38"
    somatic_dir="$somatic_dir_38"
    paths=(
      "${ref_dir}/Homo_sapiens_assembly38.fasta"
      "${ref_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz"
      "${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
      "${ref_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
      "${ref_dir}/hapmap_3.3.hg38.vcf.gz"
      "${ref_dir}/1000G_omni2.5.hg38.vcf.gz"
      "${ref_dir}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
      "${somatic_dir}/af-only-gnomad.hg38.vcf.gz"
      "${somatic_dir}/1000g_pon.hg38.vcf.gz"
    )
  fi

  echo "  Reference files (hg${build})"
  for path in "${paths[@]}"; do
    local status="MISSING"
    [ -s "$path" ] && status="OK"
    printf "    %-40s %s\n" "$(basename "$path")" "$status"
  done
  echo ""
}

refs_complete_for_build() {
  local build=$1
  local ref_dir=""
  local somatic_dir=""
  local -a paths=()

  if [ "$build" = "19" ]; then
    ref_dir="$refs_dir_19"
    somatic_dir="$somatic_dir_19"
    paths=(
      "${ref_dir}/Homo_sapiens_assembly19.fasta"
      "${ref_dir}/Homo_sapiens_assembly19.known_indels_20120518.vcf"
      "${ref_dir}/dbsnp_138.b37.vcf.gz"
      "${ref_dir}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
      "${ref_dir}/hapmap_3.3.b37.vcf.gz"
      "${ref_dir}/1000G_omni2.5.b37.vcf.gz"
      "${ref_dir}/1000G_phase1.snps.high_confidence.b37.vcf.gz"
      "${somatic_dir}/af-only-gnomad.raw.sites.vcf.gz"
      "${somatic_dir}/Mutect2-WGS-panel-b37.vcf.gz"
    )
  else
    ref_dir="$refs_dir_38"
    somatic_dir="$somatic_dir_38"
    paths=(
      "${ref_dir}/Homo_sapiens_assembly38.fasta"
      "${ref_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz"
      "${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
      "${ref_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
      "${ref_dir}/hapmap_3.3.hg38.vcf.gz"
      "${ref_dir}/1000G_omni2.5.hg38.vcf.gz"
      "${ref_dir}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
      "${somatic_dir}/af-only-gnomad.hg38.vcf.gz"
      "${somatic_dir}/1000g_pon.hg38.vcf.gz"
    )
  fi

  for path in "${paths[@]}"; do
    [ -s "$path" ] || return 1
  done
  return 0
}

run_validate_only() {
  local errors=0
  local status=""

  log "Running validate-only..."
  echo "  Mode: tools_dir=${PREFIX} refs_dir=${REFS_DIR} build=${BUILD}"
  echo ""

  if $VALIDATE_TOOLS; then
    echo "  Tools"

    status="MISSING"; [ -x "$fastqc" ] && status="OK"; [ "$status" = "MISSING" ] && errors=$((errors + 1)); printf "    %-30s %s\n" "FastQC" "$status"
    status="MISSING"; [ -f "${trim_path}/trimmomatic-0.39.jar" ] && status="OK"; [ "$status" = "MISSING" ] && errors=$((errors + 1)); printf "    %-30s %s\n" "Trimmomatic" "$status"
    status="MISSING"; [ -x "$bwamem" ] && status="OK"; [ "$status" = "MISSING" ] && errors=$((errors + 1)); printf "    %-30s %s\n" "BWA-MEM2" "$status"
    status="MISSING"; [ -x "$samtools" ] && status="OK"; [ "$status" = "MISSING" ] && errors=$((errors + 1)); printf "    %-30s %s\n" "SAMtools" "$status"
    status="MISSING"; [ -x "$gatk" ] && status="OK"; [ "$status" = "MISSING" ] && errors=$((errors + 1)); printf "    %-30s %s\n" "GATK" "$status"
    status="MISSING"; command -v bcftools &>/dev/null && status="OK"; [ "$status" = "MISSING" ] && errors=$((errors + 1)); printf "    %-30s %s\n" "bcftools" "$status"
    status="MISSING"; [ -f "${mosdepth_path}/scripts/plot-dist.py" ] && status="OK"; [ "$status" = "MISSING" ] && errors=$((errors + 1)); printf "    %-30s %s\n" "mosdepth scripts" "$status"
    status="MISSING"; [ -f "${annovar_path}/table_annovar.pl" ] && status="OK"; [ "$status" = "MISSING" ] && errors=$((errors + 1)); printf "    %-30s %s\n" "ANNOVAR" "$status"
    echo ""
  fi

  if $VALIDATE_REFS; then
    echo "  References"
    for build in "${SELECTED_BUILDS[@]}"; do
      if refs_complete_for_build "$build"; then
        printf "    %-30s %s\n" "hg${build}" "OK"
      else
        printf "    %-30s %s\n" "hg${build}" "MISSING"
        errors=$((errors + 1))
      fi
    done
    echo ""
    for build in "${SELECTED_BUILDS[@]}"; do
      print_ref_summary_for_build "$build"
    done
  fi

  if [ "$errors" -gt 0 ]; then
    echo -e "${FAIL}  validation failed: ${errors} missing item(s)"
    exit 1
  fi
  echo -e "${PASS}  validation passed"
  exit 0
}

if $VALIDATE_ONLY; then
  run_validate_only
fi

log "Checking system dependencies..."
need wget "wget"
need tar "tar"
if $DO_SOFTWARE; then
  need java "JDK >= 17"
  need conda "Miniconda/Anaconda"
  need git "git"
  need unzip "unzip"
fi
if $DO_CONDA; then
  need conda "Miniconda/Anaconda"
  need git "git"
fi
if $DO_REFS; then
  need gsutil "Google Cloud SDK"
fi
log "System dependencies OK"

ensure_install_roots

if $DO_SOFTWARE; then
  need_parallel
  log "Checking software install candidates in ${PREFIX}..."
  declare -a SW_TASK_IDS=()
  declare -a SW_TASK_LABELS=()
  queue_sw_task() {
    SW_TASK_IDS+=("$1")
    SW_TASK_LABELS+=("$2")
  }

  if [ -x "$fastqc" ]; then
    echo -e "${SKIP}  FastQC already at $fastqc"
  elif ask_permission "FastQC"; then
    queue_sw_task "fastqc" "FastQC"
  else
    echo -e "${SKIP}  FastQC (user declined)"
  fi

  if [ -f "${trim_path}/trimmomatic-0.39.jar" ]; then
    echo -e "${SKIP}  Trimmomatic already at $trim_path"
  elif ask_permission "Trimmomatic"; then
    queue_sw_task "trimmomatic" "Trimmomatic"
  else
    echo -e "${SKIP}  Trimmomatic (user declined)"
  fi

  if [ -x "$bwamem" ]; then
    echo -e "${SKIP}  BWA-MEM2 already at $bwamem"
  elif ask_permission "BWA-MEM2"; then
    queue_sw_task "bwamem2" "BWA-MEM2"
  else
    echo -e "${SKIP}  BWA-MEM2 (user declined)"
  fi

  if [ -x "$samtools" ]; then
    echo -e "${SKIP}  SAMtools already at $samtools"
  elif ask_permission "SAMtools"; then
    queue_sw_task "samtools" "SAMtools"
  else
    echo -e "${SKIP}  SAMtools (user declined)"
  fi

  if command -v bcftools &>/dev/null; then
    echo -e "${SKIP}  bcftools already in PATH ($(bcftools --version | head -1))"
  elif ask_permission "bcftools (conda base)"; then
    queue_sw_task "bcftools" "bcftools"
  else
    echo -e "${SKIP}  bcftools (user declined)"
  fi

  if [ -x "$gatk" ]; then
    echo -e "${SKIP}  GATK already at $gatk"
  elif ask_permission "GATK"; then
    queue_sw_task "gatk" "GATK"
  else
    echo -e "${SKIP}  GATK (user declined)"
  fi

  plot_script="${mosdepth_path}/scripts/plot-dist.py"
  if [ -f "$plot_script" ]; then
    echo -e "${SKIP}  mosdepth scripts already at ${mosdepth_path}/scripts/"
  elif ask_permission "mosdepth plot-dist.py"; then
    queue_sw_task "mosdepth_plot" "mosdepth plot-dist.py"
  else
    echo -e "${SKIP}  mosdepth plot-dist.py (user declined)"
  fi

  if [ "${#SW_TASK_IDS[@]}" -gt 0 ]; then
    log "Installing approved software in GNU parallel..."
    need gcc "gcc"
    need make "make"

    sw_jobs_file=$(mktemp)
    for i in "${!SW_TASK_IDS[@]}"; do
      task_id="${SW_TASK_IDS[$i]}"
      task_label="${SW_TASK_LABELS[$i]}"
      case "$task_id" in
        fastqc)
          printf "%s\t%s\n" "$task_label" "fastqc_zip='${PREFIX}/FastQC.zip'; download 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip' \"\$fastqc_zip\" && unzip -q \"\$fastqc_zip\" -d '${PREFIX}' && chmod +x '${fastqc}' && rm -f \"\$fastqc_zip\" && [ -x '${fastqc}' ]" >> "$sw_jobs_file"
          ;;
        trimmomatic)
          printf "%s\t%s\n" "$task_label" "trim_zip='${PREFIX}/Trimmomatic-0.39.zip'; download 'https://github.com/usadellab/Trimmomatic/releases/download/v0.39/Trimmomatic-0.39.zip' \"\$trim_zip\" && unzip -q \"\$trim_zip\" -d '${PREFIX}' && rm -f \"\$trim_zip\" && [ -f '${trim_path}/trimmomatic-0.39.jar' ]" >> "$sw_jobs_file"
          ;;
        bwamem2)
          printf "%s\t%s\n" "$task_label" "bwa_tar='${PREFIX}/bwa-mem2.tar.bz2'; download 'https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.3/bwa-mem2-2.3_x64-linux.tar.bz2' \"\$bwa_tar\" && tar -xjf \"\$bwa_tar\" -C '${PREFIX}' && rm -f \"\$bwa_tar\" && [ -x '${bwamem}' ]" >> "$sw_jobs_file"
          ;;
        samtools)
          printf "%s\t%s\n" "$task_label" "samtools_tar='${PREFIX}/samtools-1.22.1.tar.bz2'; samtools_src='${PREFIX}/samtools-1.22.1'; download 'https://github.com/samtools/samtools/releases/download/1.22.1/samtools-1.22.1.tar.bz2' \"\$samtools_tar\" && tar -xjf \"\$samtools_tar\" -C '${PREFIX}' && (cd \"\$samtools_src\" && ./configure --without-curses --prefix=\"\$samtools_src\" && make -j\"$(nproc)\" && make install) && rm -f \"\$samtools_tar\" && [ -x '${samtools}' ]" >> "$sw_jobs_file"
          ;;
        bcftools)
          printf "%s\t%s\n" "$task_label" "conda install -y -n base -c bioconda bcftools 2>/dev/null" >> "$sw_jobs_file"
          ;;
        gatk)
          printf "%s\t%s\n" "$task_label" "gatk_zip='${PREFIX}/gatk-4.6.2.0.zip'; download 'https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip' \"\$gatk_zip\" && unzip -q \"\$gatk_zip\" -d '${PREFIX}' && rm -f \"\$gatk_zip\" && [ -x '${gatk}' ]" >> "$sw_jobs_file"
          ;;
        mosdepth_plot)
          printf "%s\t%s\n" "$task_label" "mkdir -p '${mosdepth_path}/scripts' && download 'https://raw.githubusercontent.com/brentp/mosdepth/v0.3.3/scripts/plot-dist.py' '${plot_script}' && [ -f '${plot_script}' ]" >> "$sw_jobs_file"
          ;;
      esac
    done

    if ! parallel_eval_tab_jobs "$sw_jobs_file" 0; then
      log "Software tasks failed"
    fi
    rm -f "$sw_jobs_file"
  fi

  if [ -f "${annovar_path}/table_annovar.pl" ]; then
    echo -e "${SKIP}  ANNOVAR already at $annovar_path"
  else
    echo -e "\e[33m[MANUAL]\e[m  ANNOVAR requires free academic registration."
    echo "  1. Register at: https://annovar.openbioinformatics.org/en/latest/"
    echo "  2. Download annovar.latest.tar.gz"
    echo "  3. Extract to: $annovar_path"
    echo "  Then re-run this script."
  fi
fi

if $DO_CONDA; then
  log "Setting up conda environments..."

  create_env() {
    local env_name=$1
    shift
    if conda_env_exists "$env_name"; then
      echo -e "${SKIP}  conda env '$env_name' already exists"
    else
      log "Creating conda env: $env_name"
      conda create -y -n "$env_name" -c bioconda -c conda-forge "$@" \
        && echo -e "${PASS}  conda env '$env_name'" \
        || echo -e "${FAIL}  conda env '$env_name'"
    fi
  }

  if conda_env_exists wgs_parallel; then
    echo -e "${SKIP}  conda env 'wgs_parallel' already exists"
  else
    log "Creating conda env from ${REPO_DIR}/wgs_parallel.yml"
    conda env create -f "${REPO_DIR}/wgs_parallel.yml" \
      && echo -e "${PASS}  conda env 'wgs_parallel'" \
      || echo -e "${FAIL}  conda env 'wgs_parallel'"
  fi
fi

if $DO_REFS; then
  need_parallel
  declare -a REF_TASK_BUILDS=()
  declare -a REF_TASK_LABELS=()

  for build in "${SELECTED_BUILDS[@]}"; do
    if refs_complete_for_build "$build"; then
      echo -e "${SKIP}  hg${build} references already present"
    elif ask_permission "reference bundle hg${build}"; then
      REF_TASK_BUILDS+=("$build")
      REF_TASK_LABELS+=("hg${build} references")
    else
      echo -e "${SKIP}  hg${build} references (user declined)"
    fi
  done

  if [ "${#REF_TASK_BUILDS[@]}" -gt 0 ]; then
    log "Downloading approved references in GNU parallel..."
    ref_jobs_file=$(mktemp)
    for i in "${!REF_TASK_BUILDS[@]}"; do
      build="${REF_TASK_BUILDS[$i]}"
      printf "%s\t%s\n" "${REF_TASK_LABELS[$i]}" "download_refs_for_build '${build}'" >> "$ref_jobs_file"
    done
    if ! parallel_eval_tab_jobs "$ref_jobs_file" 0; then
      log "Reference tasks failed"
    fi
    rm -f "$ref_jobs_file"
  fi
fi

if ! python3 -c 'import yaml' >/dev/null 2>&1; then
  echo "PyYAML is required to write config.yaml; installing it now."
  if command -v conda >/dev/null 2>&1; then
    conda install -y -n base pyyaml || { echo "Unable to install PyYAML" >&2; exit 1; }
  elif python3 -m pip --version >/dev/null 2>&1; then
    python3 -m pip install --user PyYAML || { echo "Unable to install PyYAML" >&2; exit 1; }
  else
    echo "Install PyYAML for python3 and re-run setup.sh." >&2
    exit 1
  fi
fi
write_config

log "Setup complete. Updated ${REPO_DIR}/config.yaml"
echo ""
echo "  Software"
printf "    %-30s %s\n" "FastQC" "$([ -x "$fastqc" ] && echo OK || echo MISSING)"
printf "    %-30s %s\n" "Trimmomatic" "$([ -f "${trim_path}/trimmomatic-0.39.jar" ] && echo OK || echo MISSING)"
printf "    %-30s %s\n" "BWA-MEM2" "$([ -x "$bwamem" ] && echo OK || echo MISSING)"
printf "    %-30s %s\n" "SAMtools" "$([ -x "$samtools" ] && echo OK || echo MISSING)"
printf "    %-30s %s\n" "GATK" "$([ -x "$gatk" ] && echo OK || echo MISSING)"
printf "    %-30s %s\n" "bcftools" "$(command -v bcftools &>/dev/null && echo OK || echo MISSING)"
printf "    %-30s %s\n" "mosdepth scripts" "$([ -f "${mosdepth_path}/scripts/plot-dist.py" ] && echo OK || echo MISSING)"
printf "    %-30s %s\n" "ANNOVAR" "$([ -f "${annovar_path}/table_annovar.pl" ] && echo OK || echo MANUAL)"
echo ""
echo "  Conda environments"
for env in wgs_parallel; do
  printf "    %-30s %s\n" "$env" "$(conda_env_exists "$env" && echo OK || echo MISSING)"
done
echo ""
for build in "${SELECTED_BUILDS[@]}"; do
  print_ref_summary_for_build "$build"
done
