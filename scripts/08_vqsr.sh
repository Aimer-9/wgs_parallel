#!/bin/bash
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
source "${REPO_DIR}/scripts/00_util.sh"

work_dir=$1
hsa_version=$2
gatk=$3
ref_fa=$4
known_indels=$5
dbsnp=$6
Mills=$7
hapmap=$8
file_1000G_omni=$9
file_1000G_phase1=${10}
MAX_VQSR_JOBS=${11:-10}
MAX_GENOTYPE_JOBS=${12:-10}
GENOTYPE_EXTRA_ARGS=${13:-}
VQSR_SNP_TRUTH_SENSITIVITY=${14:-99.5}
VQSR_INDEL_TRUTH_SENSITIVITY=${15:-99.0}
VQSR_EXTRA_ARGS=${16:-}
bcftools=${17:-bcftools}

num_samples=$(wc -l < "${work_dir}/sample_list.txt")
[ "$MAX_VQSR_JOBS" -gt "$num_samples" ] && MAX_VQSR_JOBS=$num_samples
[ "$MAX_GENOTYPE_JOBS" -lt "$MAX_VQSR_JOBS" ] && MAX_VQSR_JOBS=$MAX_GENOTYPE_JOBS

haplotype_dir="${work_dir}/output/07_HaplotypeCaller"
vqsr_dir="${work_dir}/output/08_VQSR"
if [ ! -d "$vqsr_dir" ]; then mkdir -p "$vqsr_dir"; fi

# Full VQSR pipeline for one sample:
#   GenotypeGVCFs → VariantRecalibrator (SNP + INDEL)
#   → ApplyVQSR (SNP + INDEL) → bcftools PASS filter
run_vqsr() {
  local sample_id=$1
  local gvcf="${haplotype_dir}/${sample_id}.g.vcf.gz"
  local vcf="${vqsr_dir}/${sample_id}.vcf"
  local recal_snp="${vqsr_dir}/${sample_id}_vqsr_SNP.vcf"
  local recal_indel="${vqsr_dir}/${sample_id}_vqsr_indel.vcf"
  local applied_snp="${vqsr_dir}/${sample_id}.SNP.recalibrated.vcf"
  local applied_both="${vqsr_dir}/${sample_id}.indel.SNP.recalibrated.vcf"
  local final_vcf="${vqsr_dir}/${sample_id}.indel.SNP.recalibrated.PASS.vcf"
  local -a genotype_extra=() vqsr_extra=()
  read_extra_args "$GENOTYPE_EXTRA_ARGS" genotype_extra || return 1
  read_extra_args "$VQSR_EXTRA_ARGS" vqsr_extra || return 1

  if [ ! -s "$gvcf" ]; then
    echo "Error: ${gvcf} not found — HaplotypeCaller step may have failed"; return 1
  fi

  # GenotypeGVCFs
  if [ ! -s "$vcf" ]; then
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mGenotypeGVCFs: ${sample_id} processing\e[m"
    "$gatk" GenotypeGVCFs \
      --verbosity WARNING \
      -R "$ref_fa" \
      --variant "$gvcf" \
      "${genotype_extra[@]}" \
      -O "$vcf"
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mGenotypeGVCFs: ${sample_id} done\e[m"
  else
    echo -e "\e[37;42mGenotypeGVCFs: ${sample_id} exists\e[m"
  fi

  if [ ! -s "$vcf" ]; then
    echo "Error: genotyped VCF missing for ${sample_id}"; return 1
  fi

  # VariantRecalibrator — SNP
  if [ ! -s "$recal_snp" ]; then
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mVQSR SNP recal: ${sample_id} processing\e[m"
    "$gatk" VariantRecalibrator \
      --verbosity WARNING \
      -R "$ref_fa" \
      -V "$vcf" \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0   "$hapmap" \
      --resource:omni,known=false,training=true,truth=false,prior=12.0    "$file_1000G_omni" \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0   "$file_1000G_phase1" \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0    "$dbsnp" \
      -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
      "${vqsr_extra[@]}" \
      -mode SNP \
      -O "$recal_snp" \
      --tranches-file "${vqsr_dir}/${sample_id}_output_SNP.tranches" \
      --rscript-file  "${vqsr_dir}/${sample_id}_output_SNP.plots.R"
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mVQSR SNP recal: ${sample_id} done\e[m"
  else
    echo -e "\e[37;42mVQSR SNP recal: ${sample_id} exists\e[m"
  fi

  # VariantRecalibrator — INDEL
  if [ ! -s "$recal_indel" ]; then
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mVQSR INDEL recal: ${sample_id} processing\e[m"
    "$gatk" VariantRecalibrator \
      --verbosity WARNING \
      -R "$ref_fa" \
      -V "$vcf" \
      --resource:mills,known=false,training=true,truth=true,prior=12.0  "$Mills" \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0  "$dbsnp" \
      -an QD -an DP -an FS -an ReadPosRankSum -an MQRankSum \
      "${vqsr_extra[@]}" \
      -mode INDEL \
      -O "$recal_indel" \
      --tranches-file "${vqsr_dir}/${sample_id}_output_indel.tranches" \
      --rscript-file  "${vqsr_dir}/${sample_id}_output_indel.plots.R"
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mVQSR INDEL recal: ${sample_id} done\e[m"
  else
    echo -e "\e[37;42mVQSR INDEL recal: ${sample_id} exists\e[m"
  fi

  if [ ! -s "$recal_snp" ] || [ ! -s "$recal_indel" ]; then
    echo "Error: recalibration file(s) missing for ${sample_id}"; return 1
  fi

  # ApplyVQSR — SNP
  if [ ! -s "$applied_snp" ]; then
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mApplyVQSR SNP: ${sample_id} processing\e[m"
    "$gatk" ApplyVQSR \
      --verbosity WARNING \
      -V "$vcf" \
      --recal-file "$recal_snp" \
      -mode SNP \
      --tranches-file "${vqsr_dir}/${sample_id}_output_SNP.tranches" \
      --truth-sensitivity-filter-level "$VQSR_SNP_TRUTH_SENSITIVITY" \
      -O "$applied_snp"
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mApplyVQSR SNP: ${sample_id} done\e[m"
  else
    echo -e "\e[37;42mApplyVQSR SNP: ${sample_id} exists\e[m"
  fi

  # ApplyVQSR — INDEL (applied on top of SNP-recalibrated VCF)
  if [ ! -s "$applied_both" ]; then
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mApplyVQSR INDEL: ${sample_id} processing\e[m"
    "$gatk" ApplyVQSR \
      --verbosity WARNING \
      -V "$applied_snp" \
      --recal-file "$recal_indel" \
      -mode INDEL \
      --tranches-file "${vqsr_dir}/${sample_id}_output_indel.tranches" \
      --truth-sensitivity-filter-level "$VQSR_INDEL_TRUTH_SENSITIVITY" \
      -O "$applied_both"
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mApplyVQSR INDEL: ${sample_id} done\e[m"
  else
    echo -e "\e[37;42mApplyVQSR INDEL: ${sample_id} exists\e[m"
  fi

  if [ ! -s "$applied_both" ]; then
    echo "Error: VQSR-applied VCF missing for ${sample_id}"; return 1
  fi

  # Filter PASS variants only
  if [ ! -s "$final_vcf" ]; then
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mPASS filter: ${sample_id} processing\e[m"
    "$bcftools" view -f PASS "$applied_both" -o "$final_vcf"
    date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mPASS filter: ${sample_id} done\e[m"
  else
    echo -e "\e[37;42mPASS filter: ${sample_id} exists\e[m"
  fi
}

cleanup_run_vqsr() {
  local sample_id=$1
  rm -f "${vqsr_dir}/${sample_id}.vcf"* \
    "${vqsr_dir}/${sample_id}_vqsr_SNP.vcf"* \
    "${vqsr_dir}/${sample_id}_vqsr_indel.vcf"* \
    "${vqsr_dir}/${sample_id}.SNP.recalibrated.vcf"* \
    "${vqsr_dir}/${sample_id}.indel.SNP.recalibrated.vcf"* \
    "${vqsr_dir}/${sample_id}.indel.SNP.recalibrated.PASS.vcf"* \
    "${vqsr_dir}/${sample_id}_output_SNP.tranches" \
    "${vqsr_dir}/${sample_id}_output_SNP.plots.R" \
    "${vqsr_dir}/${sample_id}_output_indel.tranches" \
    "${vqsr_dir}/${sample_id}_output_indel.plots.R"
}

date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mVQSR: start\e[m"
parallel_run_sample_list "${work_dir}/sample_list.txt" "$MAX_VQSR_JOBS" run_vqsr vqsr || exit 1
date +"%Y-%m-%d %H:%M:%S" && echo -e "\e[37;42mVQSR: all done\e[m"
