#!/bin/bash
set -e
tissue=$1
num=$2

  osca="/soft/osca-0.46.1-linux-x86_64/osca"
  res_dir="/juice_data/zhh_independent/ZHH_update/vQTL/res/${tissue}/zscore_update"
  genofile="/juice_data/zhh_independent/ZHH_update/Genotype/${tissue}/${tissue}.autosomal.geno"
  phenofile="/juice_data/zhh_independent/ZHH_update/vQTL/pheno/${tissue}/${tissue}.expr_tmm_inv.zscore"

mkdir -p ${res_dir}

#step2 vQTL analysis
echo "Tasks2 pocessing..."
  if [[ -s "${res_dir}/${tissue}.vqtl_chr${num}.all.query.gz" ]];then exit;fi

  $osca --vqtl --vqtl-mtd 2 \
    --bfile ${genofile} \
    --befile ${phenofile} \
    --cis-wind 1000 \
    --chr ${num} \
    --thread-num 10 \
    --out ${res_dir}/${tissue}.vqtl_chr${num}

  #$osca --beqtl-summary "${res_dir}/${tissue}.vqtl_chr${num}" --query 1.0e-5 --out "${res_dir}/${tissue}.vqtl_chr${num}.query"
  $osca --beqtl-summary "${res_dir}/${tissue}.vqtl_chr${num}" --query 1 --out "${res_dir}/${tissue}.vqtl_chr${num}.all.query"
 
    gzip "${res_dir}/${tissue}.vqtl_chr${num}.all.query"
    gzip "${res_dir}/${tissue}.vqtl_chr${num}.esi"
    #rm "${res_dir}/${tissue}.vqtl_chr${num}.query_1_1.log"
    rm "${res_dir}/${tissue}.vqtl_chr${num}.all.query_1_1.log"

  #  echo "Processed and deleted: $file"
  #else
  #  echo "Warning: ${res_dir}/${tissue}.vqtl_chr${num}.query not generated, skipping."
  #fi

echo "Tasks2 completed."
