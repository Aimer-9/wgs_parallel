#!/bin/bash
# Stages 09 and 11: ANNOVAR conversion, annotation, and summary statistics.
# The same script annotates germline VQSR output and somatic Mutect2 output;
# wgs.sh supplies a different VCF directory, suffix, sample list, and output dir.
#
# Positional arguments:
#   1 work_dir; 2 reference build; 3 ANNOVAR installation;
#   4 input VCF directory; 5 per-sample VCF suffix; 6 CPU fraction;
#   7 update ClinVar (yes/no); 8 maximum concurrent samples;
#   9 sample-list file; 10 protocol template containing {clinvar};
#  11 ANNOVAR operation list; 12 extra table_annovar arguments;
#  13 annotation output directory; 14 Conda executable.
#
# Required ANNOVAR databases are checked before parallel sample work. ClinVar is
# resolved dynamically to the newest local version, optionally refreshed using
# update_annovar_db. Per sample, VCF is converted to avinput and table_annovar is
# run with calculated threads. Interrupted conversion/annotation files are
# cleaned through checkpoints. Final tables include functional, exonic, and
# ClinVar-focused summaries.
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "${REPO_DIR}/scripts/00_util.sh"

work_dir=$1
hsa_version=$2
annovar_path=$3
vcf_prefix=$4        # directory containing per-sample VCF files
vcf_suffix=$5        # filename suffix, e.g. ".indel.SNP.recalibrated.PASS.vcf"
MAX_PROCESSOR_USE_PERCENT=$6
update_clinvar=$7    # "yes" or "no"
MAX_ANNOVAR_JOBS=${8:-10}
sample_list_file=${9:-${work_dir}/sample_list.txt}
ANNOTATION_PROTOCOL_TEMPLATE=${10:-refGene,cytoBand,avsnp150,exac03,{clinvar}}
OPERATION_TYPE=${11:-g,r,f,f,f}
ANNOVAR_EXTRA_ARGS=${12:-}
annovar_dir=${13:-${vcf_prefix}/annovar}
conda_bin=${14:-conda}

# Older generated configs could append one extra closing brace while passing
# the placeholder through positional arguments. Normalize only that known
# typo; all other brace patterns remain hard errors below.
ANNOTATION_PROTOCOL_TEMPLATE=${ANNOTATION_PROTOCOL_TEMPLATE//\{clinvar\}\}/\{clinvar\}}

# ANNOVAR requires one operation code for every protocol entry. Validate again
# here so direct script calls receive the same protection as config.py users.
IFS=',' read -r -a protocol_items <<< "$ANNOTATION_PROTOCOL_TEMPLATE"
IFS=',' read -r -a operation_items <<< "$OPERATION_TYPE"
for protocol_item in "${protocol_items[@]}"; do
  if { [[ "$protocol_item" == *'{'* ]] || [[ "$protocol_item" == *'}'* ]]; } \
    && [ "$protocol_item" != '{clinvar}' ]; then
    echo "Error: malformed ANNOVAR protocol placeholder: $protocol_item" >&2
    echo "Only the exact token {clinvar} may contain braces." >&2
    exit 1
  fi
done
if [ "${#protocol_items[@]}" -ne "${#operation_items[@]}" ]; then
  echo "Error: ANNOVAR protocol/operation count mismatch: ${#protocol_items[@]} protocols but ${#operation_items[@]} operations" >&2
  echo "  protocol : $ANNOTATION_PROTOCOL_TEMPLATE" >&2
  echo "  operation: $OPERATION_TYPE" >&2
  exit 1
fi

if [ ! -s "$sample_list_file" ]; then
  echo "No samples available for annotation: $sample_list_file"
  exit 0
fi
num_samples=$(wc -l < "$sample_list_file")
nproc=$(nproc)
[ "$MAX_ANNOVAR_JOBS" -gt "$num_samples" ] && MAX_ANNOVAR_JOBS=$num_samples
threads_per_sample=$(calculate_threads "$MAX_PROCESSOR_USE_PERCENT" "$MAX_ANNOVAR_JOBS")

humandb="${annovar_path}/humandb"
if [ ! -d "$humandb" ]; then mkdir -p "$humandb"; fi
if [ ! -d "$annovar_dir" ]; then mkdir -p "$annovar_dir"; fi

# Ensure stable protocol resources exist; missing databases are downloaded into
# ANNOVAR humandb using the selected hg19/hg38 build prefix.
ensure_db() {
  local db_name=$1
  local db_file="${humandb}/hg${hsa_version}_${db_name}.txt"
  if [ ! -s "$db_file" ]; then
    date +"%Y-%m-%d %H:%M:%S" && echo "Downloading ANNOVAR database: ${db_name}"
    perl "${annovar_path}/annotate_variation.pl" \
      -buildver "hg${hsa_version}" -downdb -webfrom annovar \
      "$db_name" "$humandb"
  fi
}

ensure_db refGene
ensure_db cytoBand
ensure_db exac03
ensure_db avsnp150

# Map the numeric pipeline build to update_annovar_db's GRCh naming convention.
if [ "${hsa_version}" -eq 19 ]; then
  clinvar_build=GRCh37
elif [ "${hsa_version}" -eq 38 ]; then
  clinvar_build=GRCh38
fi

update_resources_py="${REPO_DIR}/scripts/update_resources.py"
avinput_converter_py="${REPO_DIR}/scripts/avinput2annovardb.py"

if [ "$update_clinvar" = "yes" ]; then
  if [ ! -f "$update_resources_py" ] || [ ! -f "$avinput_converter_py" ]; then
    echo "Error: integrated ClinVar updater scripts are missing under ${REPO_DIR}/scripts" >&2
    exit 1
  fi
  date +"%Y-%m-%d %H:%M:%S" && echo "Updating ClinVar database (${clinvar_build})..."
  "$conda_bin" run --no-capture-output -n wgs_parallel \
    python "$update_resources_py" \
    -d clinvar \
    -hp "$humandb" \
    -a  "$annovar_path" \
    -g  "$clinvar_build"
fi

# Select the most recently installed ClinVar file. Restrict the release suffix
# to eight digits so stale files such as hg19_clinvar_20250120}.txt cannot be
# expanded into an invalid ANNOVAR protocol name.
latest_clinvar_file=""
latest_clinvar_mtime=0
while IFS= read -r candidate; do
  candidate_name=$(basename "$candidate")
  candidate_version=${candidate_name#hg${hsa_version}_clinvar_}
  candidate_version=${candidate_version%.txt}
  if [[ "$candidate_version" =~ ^[0-9]{8}$ ]]; then
    candidate_mtime=$(stat -c '%Y' "$candidate" 2>/dev/null || echo 0)
    if [ "$candidate_mtime" -ge "$latest_clinvar_mtime" ]; then
      latest_clinvar_file="$candidate"
      latest_clinvar_mtime=$candidate_mtime
    fi
  fi
done < <(find "$humandb" -maxdepth 1 -type f -name "hg${hsa_version}_clinvar_*.txt" -print)

if [ -z "$latest_clinvar_file" ]; then
  echo "Error: no ClinVar database found in ${humandb} for hg${hsa_version}"
  echo "Re-run with update_clinvar=\"yes\" to download it."
  exit 1
fi

clinvar_version=$(basename -s .txt "$latest_clinvar_file" | sed "s/hg${hsa_version}_//")
date +"%Y-%m-%d %H:%M:%S" && echo "Using ClinVar version: ${clinvar_version}"

# Convert each input VCF to ANNOVAR's avinput representation. Mutect2 inputs use
# multi-sample conversion so all tumor-normal loci remain in one avinput file.
convert_to_avinput() {
  local sample_id=$1
  local vcf="${vcf_prefix}/${sample_id}${vcf_suffix}"
  local avinput="${annovar_dir}/${sample_id}.avinput"
  local multisample_marker="${avinput}.multisample"
  local -a conversion_args=()

  if [[ "$vcf_suffix" == *mutect2*.vcf* ]]; then
    conversion_args+=("-allsample" "-withfreq")
  fi

  if [ ! -s "$vcf" ]; then
    echo "Warning: VCF not found for ${sample_id} — skipping (${vcf})"; return 0
  fi
  if [ -s "$avinput" ]; then
    if [[ "$vcf_suffix" != *mutect2*.vcf* ]] || [ -e "$multisample_marker" ]; then
      echo -e "\e[37;42mannovar convert: ${sample_id} exists\e[m"; return 0
    fi
    echo "Legacy somatic avinput detected for ${sample_id}; rebuilding"
    rm -f "$avinput" "$multisample_marker" \
      "${annovar_dir}/${sample_id}.hg${hsa_version}_"* \
      "${annovar_dir}/${sample_id}.log" \
      "${annovar_dir}/${sample_id}.invalid_input"
  fi

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar convert: ${sample_id} processing\e[m"
  perl "${annovar_path}/convert2annovar.pl" \
    -format vcf4 \
    -includeinfo \
    "${conversion_args[@]}" \
    "$vcf" \
    -out "$avinput" || return 1
  [ -s "$avinput" ] || return 1
  if [[ "$vcf_suffix" == *mutect2*.vcf* ]]; then
    : > "$multisample_marker"
  fi
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar convert: ${sample_id} done\e[m"
}

cleanup_convert_to_avinput() {
  rm -f "${annovar_dir}/$1.avinput" "${annovar_dir}/$1.avinput.multisample"
}

date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar convert: start\e[m"
parallel_run_sample_list "$sample_list_file" "$MAX_ANNOVAR_JOBS" convert_to_avinput annovar_convert || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar convert: all done\e[m"

# Replace {clinvar} only after the latest installed database version is known.
# Protocol entries and operation types stay aligned by their comma positions.
ANNOTATION_PROTOCOL=${ANNOTATION_PROTOCOL_TEMPLATE//\{clinvar\}/$clinvar_version}
if [[ "$ANNOTATION_PROTOCOL" == *'{'* ]] || [[ "$ANNOTATION_PROTOCOL" == *'}'* ]]; then
  echo "Error: unresolved brace in ANNOVAR protocol after ClinVar expansion: $ANNOTATION_PROTOCOL" >&2
  exit 1
fi

run_annotation() {
  local sample_id=$1
  local avinput="${annovar_dir}/${sample_id}.avinput"
  local vcf="${vcf_prefix}/${sample_id}${vcf_suffix}"
  local avoutput_prefix="${annovar_dir}/${sample_id}"
  local multianno="${avoutput_prefix}.hg${hsa_version}_multianno.txt"
  local -a extra_args=()
  read_extra_args "$ANNOVAR_EXTRA_ARGS" extra_args || return 1

  if [ -e "${avoutput_prefix}.vcfinput" ]; then
    echo "Direct-VCF somatic annotation detected for ${sample_id}; rebuilding from avinput"
    rm -f "${avoutput_prefix}.hg${hsa_version}_"* \
      "${avoutput_prefix}.log" "${avoutput_prefix}.invalid_input" \
      "${avoutput_prefix}.vcfinput"
  fi

  if [ ! -s "$avinput" ]; then
    echo "Warning: avinput not found for ${sample_id} — skipping"; return 0
  fi
  if [ -s "$multianno" ]; then
    python3 "${REPO_DIR}/scripts/rename_annovar_otherinfo.py" \
      --multianno "$multianno" --vcf "$vcf" --avinput "$avinput" || return 1
    echo -e "\e[37;42mannovar annotate: ${sample_id} exists\e[m"; return 0
  fi

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar annotate: ${sample_id} processing\e[m"
  perl "${annovar_path}/table_annovar.pl" \
    "$avinput" \
    "$humandb" \
    -buildver "hg${hsa_version}" \
    --otherinfo \
    --maxgenethread "$threads_per_sample" \
    --thread        "$threads_per_sample" \
    -out            "$avoutput_prefix" \
    -protocol       "$ANNOTATION_PROTOCOL" \
    -operation      "$OPERATION_TYPE" \
    -nastring NA \
    "${extra_args[@]}" \
    -remove || return 1
  [ -s "$multianno" ] || return 1
  python3 "${REPO_DIR}/scripts/rename_annovar_otherinfo.py" \
    --multianno "$multianno" --vcf "$vcf" --avinput "$avinput" || return 1
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar annotate: ${sample_id} done\e[m"
}

cleanup_run_annotation() {
  local prefix="${annovar_dir}/$1"
  rm -f "${prefix}.hg${hsa_version}_"* "${prefix}.log" "${prefix}.invalid_input"
}

date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar annotate: start\e[m"
parallel_run_sample_list "$sample_list_file" "$MAX_ANNOVAR_JOBS" run_annotation annovar_annotate || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar annotate: all done\e[m"

# Generate lightweight summaries from multianno columns after every parallel
# annotation is complete. These files are independently resumable.
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar stats: start\e[m"
while IFS= read -r sample_id; do
  local_multianno="${annovar_dir}/${sample_id}.hg${hsa_version}_multianno.txt"
  if [ ! -s "$local_multianno" ]; then
    echo "Warning: multianno file not found for ${sample_id} — skipping stats"; continue
  fi

  exonic_stats="${annovar_dir}/${sample_id}.hg${hsa_version}_multianno.exonic_stats.txt"
  alltype_stats="${annovar_dir}/${sample_id}.hg${hsa_version}_multianno.alltype_stats.txt"
  clinvar_out="${annovar_dir}/${sample_id}.hg${hsa_version}_multianno.clinvar.txt"

  # Count coding consequences from ExonicFunc.refGene (multianno column 9).
  if [ ! -s "$exonic_stats" ]; then
    awk -F'\t' '
      BEGIN { OFS="\t"; print "ExonicFunc.refGene", "Count" }
      NR > 1 && $9 != "" && $9 != "." {
        count[$9]++
      }
      END { for (t in count) print t, count[t] }
    ' "$local_multianno" > "$exonic_stats"
  fi

  # Count broader genomic regions from Func.refGene (multianno column 6).
  if [ ! -s "$alltype_stats" ]; then
    awk -F'\t' '
      BEGIN { OFS="\t"; print "Func.refGene", "Count" }
      NR > 1 && $6 != "" && $6 != "." {
        count[$6]++
      }
      END { for (t in count) print t, count[t] }
    ' "$local_multianno" > "$alltype_stats"
  fi

  # Retain the header and rows with a non-NA CLNALLELEID (column 21).
  if [ ! -s "$clinvar_out" ]; then
    awk -F'\t' 'NR==1 || $21 != "NA"' "$local_multianno" > "$clinvar_out"
  fi

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar stats: ${sample_id} done\e[m"
done < "$sample_list_file"
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar: all done\e[m"
