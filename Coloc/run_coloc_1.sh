#!/bin/bash
start=$1
end=$2
trait=$3
locidr="/ossfs/GWAS/${trait}/colocation/"

source /soft/miniconda3/etc/profile.d/conda.sh
conda activate myenv

   echo " processing coloc analysis between vQTL and GWAS"
script="/soft/zhh/vQTL/Coloc/coloc.GWAS_QTL_1.R"
log="${trait}.coloc_vQTL.${start}.log"

Rscript ${script} "$start" "$end" "${trait}" "${locidr}" > ./logs/${log} 2>&1

echo "Coloc analysis started, check ${log} for details."
