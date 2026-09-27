#!/bin/bash
set -e
set -x
trait=$1

locidr="/ossfs/GWAS/${trait}/colocation/"
script="/soft/zhh/vQTL/Coloc/coloc.GWAS_QTL_1.R"
source /soft/miniconda3/etc/profile.d/conda.sh
conda activate myenv

max_jobs=2
current_jobs=0

for i in $(seq 13 32); do
log="${trait}.coloc_vQTL.${i}.log"

echo " processing coloc analysis on $trait task: $i"
nohup Rscript ${script} "$i" "$i" "${trait}" "${locidr}" > ./logs/${log} 2>&1 &

    ((current_jobs += 1))

    if ((current_jobs >= max_jobs)); then
        wait -n
        ((current_jobs--))
    fi

done

wait
echo "All tasks completed."
