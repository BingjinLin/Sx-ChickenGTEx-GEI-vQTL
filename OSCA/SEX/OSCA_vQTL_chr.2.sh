#!/bin/bash
set -e
tissue=$1
SEX=$2
num=$3

  osca="/soft/osca-0.46.1-linux-x86_64/osca"
  res_dir="/juice_data/zhh_independent/ZHH_update/vQTL/res/SEX/${tissue}/${SEX}"
  genofile="/juice_data/zhh_independent/ZHH_update/Genotype/${tissue}/${tissue}_${SEX}.autosomal.geno"
  phenofile="/juice_data/zhh_independent/ZHH_update/vQTL/pheno/${tissue}/${tissue}_${SEX}.expr_tmm_inv.zscore"

mkdir -p ${res_dir}

#step2 vQTL analysis
echo "Tasks2 pocessing..."

  query_file="${res_dir}/${tissue}_${SEX}.vqtl_chr${num}.all.query"
  esi_file="${res_dir}/${tissue}_${SEX}.vqtl_chr${num}.esi"

  $osca --vqtl --vqtl-mtd 2 \
    --bfile ${genofile} \
    --befile ${phenofile} \
    --cis-wind 1000 \
    --chr ${num} \
    --thread-num 10 \
    --out ${res_dir}/${tissue}_${SEX}.vqtl_chr${num}

  $osca --beqtl-summary "${res_dir}/${tissue}_${SEX}.vqtl_chr${num}" --query 1 --out "${query_file}"


    if [[ -s "$query_file" && -s "$esi_file" ]];then

    gzip -f "${query_file}"
    rm -f "${esi_file}"
    rm -f "${res_dir}/${tissue}_${SEX}.vqtl_chr${num}.besd"
    rm -f "${res_dir}/${tissue}_${SEX}.vqtl_chr${num}.all.query_1_1.log"

    fi

echo "Tasks2 completed."
