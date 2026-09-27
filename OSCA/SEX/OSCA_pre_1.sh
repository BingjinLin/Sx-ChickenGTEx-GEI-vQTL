#!/bin/bash
set -e
tissue=$1
SEX=$2


  osca="/soft/osca-0.46.1-linux-x86_64/osca"
  res_dir="/data/vQTL/SEX/${tissue}"
  datadir="/juice_data/zhh_independent/ZHH_update/vQTL/pheno"
  
  
  phenofile="${datadir}/${tissue}/${tissue}_${SEX}.expr_tmm_inv.zscore"
  anno1="${datadir}/${tissue}/${tissue}.opi"
  anno2="${datadir}/${tissue}/${tissue}_${SEX}.oii"

#mkdir -p ${res_dir}

#step0 prepare genotype file

  genofile="/juice_data/zhh_independent/ZHH_update/Genotype/${tissue}/${tissue}.autosomal.geno"
  sample_list="/juice_data/zhh_independent/ZHH_update/Genotype/${tissue}/${tissue}_${SEX}_sample.list"
new_geno="/juice_data/zhh_independent/ZHH_update/Genotype/${tissue}/${tissue}_${SEX}.autosomal.geno"

if [ -f "$sample_list" ]; then
plink --bfile ${genofile} --keep ${sample_list} \
--recode --keep-allele-order --make-bed --const-fid \
--chr-set 39 --out ${new_geno}
else
    echo "There is no file: ${sample_list}"
fi

#step1 generate phenotype file
echo "Launching task1 for $tissue time=$SEX..."
if [ -f "${phenofile}" ]; then

mv ${phenofile}.oii ${anno2}

  $osca --efile ${phenofile} --gene-expression --make-bod --out ${phenofile}
  $osca --befile ${phenofile} --update-opi ${anno1}

mv ${anno2} ${phenofile}.oii
else
    echo "There is no file: ${phenofile}"
fi
echo "Tasks1 completed."