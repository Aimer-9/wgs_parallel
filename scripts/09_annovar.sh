#!/bin/bash
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
ANNOTATION_PROTOCOL_TEMPLATE=${10:-refGene,cytoBand,{clinvar},avsnp150,exac03}
OPERATION_TYPE=${11:-g,r,f,f,f}
ANNOVAR_EXTRA_ARGS=${12:-}
annovar_dir=${13:-${vcf_prefix}/annovar}

if [ ! -s "$sample_list_file" ]; then
  echo "No samples available for annotation: $sample_list_file"
  exit 0
fi
num_samples=$(wc -l < "$sample_list_file")
nproc=$(nproc)
[ "$MAX_ANNOVAR_JOBS" -gt "$num_samples" ] && MAX_ANNOVAR_JOBS=$num_samples
threads_per_sample=$(calculate_threads "$MAX_PROCESSOR_USE_PERCENT" "$MAX_ANNOVAR_JOBS")

humandb="${annovar_path}/humandb"
if [ ! -d "$annovar_dir" ]; then mkdir -p "$annovar_dir"; fi

# ── Database check + auto-download ───────────────────────────────────────────
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

# ── ClinVar version resolution ────────────────────────────────────────────────
# Map hsa_version integer to the genome build string used by update_annovar_db
if [ "${hsa_version}" -eq 19 ]; then
  clinvar_build=GRCh37
elif [ "${hsa_version}" -eq 38 ]; then
  clinvar_build=GRCh38
fi

update_annovar_db_dir="${annovar_path}/update_annovar_db"

# Clone update_annovar_db repo if not present
if [ ! -f "${update_annovar_db_dir}/update_resources.py" ]; then
  git clone https://github.com/mobidic/update_annovar_db.git "$update_annovar_db_dir"
  conda env create -f "${update_annovar_db_dir}/environment.yml"
fi

if [ "$update_clinvar" = "yes" ]; then
  date +"%Y-%m-%d %H:%M:%S" && echo "Updating ClinVar database (${clinvar_build})..."
  source activate update_annovar_db
  python "${update_annovar_db_dir}/update_resources.py" \
    -d clinvar \
    -hp "$humandb" \
    -a  "$annovar_path" \
    -g  "$clinvar_build"
  conda deactivate

  # Copy the downloaded file into humandb/ with the hgXX_ prefix
  latest_clinvar_txt=$(find "${update_annovar_db_dir}/clinvar/${clinvar_build}" \
    -maxdepth 1 -name 'clinvar_*.txt' -printf '%T@ %p\n' \
    | sort -n | tail -1 | awk '{print $2}')
  cp "$latest_clinvar_txt" "${humandb}/hg${hsa_version}_$(basename "$latest_clinvar_txt")"
fi

# Select the most recently downloaded local ClinVar file
latest_clinvar_file=$(find "$humandb" -maxdepth 1 \
  -name "hg${hsa_version}_clinvar_*.txt" -printf '%T@ %p\n' \
  | sort -n | tail -1 | awk '{print $2}')

if [ -z "$latest_clinvar_file" ]; then
  echo "Error: no ClinVar database found in ${humandb} for hg${hsa_version}"
  echo "Re-run with update_clinvar=\"yes\" to download it."
  exit 1
fi

clinvar_version=$(basename -s .txt "$latest_clinvar_file" | sed "s/hg${hsa_version}_//")
date +"%Y-%m-%d %H:%M:%S" && echo "Using ClinVar version: ${clinvar_version}"

# ── VCF → ANNOVAR input conversion ───────────────────────────────────────────
convert_to_avinput() {
  local sample_id=$1
  local vcf="${vcf_prefix}/${sample_id}${vcf_suffix}"
  local avinput="${annovar_dir}/${sample_id}.avinput"

  if [ ! -s "$vcf" ]; then
    echo "Warning: VCF not found for ${sample_id} — skipping (${vcf})"; return 0
  fi
  if [ -s "$avinput" ]; then
    echo -e "\e[37;42mannovar convert: ${sample_id} exists\e[m"; return 0
  fi

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar convert: ${sample_id} processing\e[m"
  perl "${annovar_path}/convert2annovar.pl" \
    -format vcf4 \
    -includeinfo \
    "$vcf" \
    -out "$avinput"
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar convert: ${sample_id} done\e[m"
}

cleanup_convert_to_avinput() {
  rm -f "${annovar_dir}/$1.avinput"
}

date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar convert: start\e[m"
parallel_run_sample_list "$sample_list_file" "$MAX_ANNOVAR_JOBS" convert_to_avinput annovar_convert || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar convert: all done\e[m"

# ── table_annovar annotation ──────────────────────────────────────────────────
# Protocol: refGene (gene), cytoBand (region), ClinVar + avsnp150 + exac03 (filter)
ANNOTATION_PROTOCOL=${ANNOTATION_PROTOCOL_TEMPLATE//\{clinvar\}/$clinvar_version}

run_annotation() {
  local sample_id=$1
  local avinput="${annovar_dir}/${sample_id}.avinput"
  local avoutput_prefix="${annovar_dir}/${sample_id}"
  local -a extra_args=()
  read_extra_args "$ANNOVAR_EXTRA_ARGS" extra_args || return 1

  if [ ! -s "$avinput" ]; then
    echo "Warning: avinput not found for ${sample_id} — skipping"; return 0
  fi
  if [ -s "${avoutput_prefix}.hg${hsa_version}_multianno.txt" ]; then
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
    -remove
  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar annotate: ${sample_id} done\e[m"
}

cleanup_run_annotation() {
  local prefix="${annovar_dir}/$1"
  rm -f "${prefix}.hg${hsa_version}_"* "${prefix}.log" "${prefix}.invalid_input"
}

date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar annotate: start\e[m"
parallel_run_sample_list "$sample_list_file" "$MAX_ANNOVAR_JOBS" run_annotation annovar_annotate || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar annotate: all done\e[m"

# ── Per-sample statistics ─────────────────────────────────────────────────────
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar stats: start\e[m"
while IFS= read -r sample_id; do
  local_multianno="${annovar_dir}/${sample_id}.hg${hsa_version}_multianno.txt"
  if [ ! -s "$local_multianno" ]; then
    echo "Warning: multianno file not found for ${sample_id} — skipping stats"; continue
  fi

  exonic_stats="${annovar_dir}/${sample_id}.hg${hsa_version}_multianno.exonic_stats.txt"
  alltype_stats="${annovar_dir}/${sample_id}.hg${hsa_version}_multianno.alltype_stats.txt"
  clinvar_out="${annovar_dir}/${sample_id}.hg${hsa_version}_multianno.clinvar.txt"

  # Exonic variant type counts (col 9 = ExonicFunc.refGene)
  if [ ! -s "$exonic_stats" ]; then
    awk -F'\t' '
      BEGIN { OFS="\t"; print "ExonicFunc.refGene", "Count" }
      NR > 1 && $9 != "" && $9 != "." {
        count[$9]++
      }
      END { for (t in count) print t, count[t] }
    ' "$local_multianno" > "$exonic_stats"
  fi

  # All functional region counts (col 6 = Func.refGene)
  if [ ! -s "$alltype_stats" ]; then
    awk -F'\t' '
      BEGIN { OFS="\t"; print "Func.refGene", "Count" }
      NR > 1 && $6 != "" && $6 != "." {
        count[$6]++
      }
      END { for (t in count) print t, count[t] }
    ' "$local_multianno" > "$alltype_stats"
  fi

  # ClinVar-annotated variants only (col 12 = CLNALLELEID, non-NA)
  if [ ! -s "$clinvar_out" ]; then
    awk -F'\t' 'NR==1 || $12 != "NA"' "$local_multianno" > "$clinvar_out"
  fi

  date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar stats: ${sample_id} done\e[m"
done < "$sample_list_file"
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mannovar: all done\e[m"
