#!/bin/bash
set -e
tissue=$1

  osca="/soft/osca-0.46.1-linux-x86_64/osca"
  res_dir="/juice_data/zhh_independent/ZHH_update/vQTL/res/${tissue}"
  genofile="/juice_data/zhh_independent/ZHH_update/Genotype/${tissue}/${tissue}.autosomal.geno"
  datadir="/juice_data/zhh_independent/ZHH_update/vQTL/pheno"
  phenofile="${datadir}/${tissue}/${tissue}.expr_tmm_inv.zscore"
  anno1="${datadir}/${tissue}/${tissue}.opi"
  anno2="${datadir}/${tissue}/${tissue}.oii"

#mkdir -p ${res_dir}

#step1 generate phenotype file
echo "Tasks1 pocessing..."

  $osca --efile ${phenofile} --gene-expression --make-bod --out ${phenofile}
  $osca --befile ${phenofile} --update-opi ${anno1}
mv ${anno2} ${phenofile}.oii

echo "Tasks1 completed."
