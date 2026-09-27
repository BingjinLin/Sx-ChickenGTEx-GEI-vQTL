#!/bin/bash
set -e
set -x

geno_file="/juice_data/zhh_independent/ZHH_update/Genotype/Chicken.SNP_Autosomal.recode"
snp_file="/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/top_pairs/multi.top_vqtl.list"
Otfile="/juice_data/zhh_independent/ZHH_update/vQTL/res/multi.top_vqtl"

#echo "Extract top vqtl"
#plink --bfile ${geno_file} --chr-set 39 --keep-allele-order --make-bed --extract ${snp_file} --out ${Otfile}
echo "Calculate LD top vqtl"
plink --bfile ${Otfile} --chr-set 39 --r2 \
  --ld-window 999999 \
  --ld-window-kb 999999 \
  --ld-window-r2 0 \
  --out ${Otfile}

