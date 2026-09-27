#!/bin/bash
set -e

source /soft/miniconda3/etc/profile.d/conda.sh
conda activate myenv
script="/soft/zhh/vQTL/Coloc/coloc_cievQTL.R"

max_jobs=4
current_jobs=0

tissue_list=($(cat /soft/zhh/vQTL/Coloc/coloc_cie_tissue.task.list))
#tissue_list=("${tissue_list[@]:10:32}")

for tis in "${tissue_list[@]}"
do
log="./coloc_log/${tis}.coloc_cieQTL.log"
echo "Coloc analysis on ${tis} started, check ${log} for details."

nohup Rscript ${script} "$tis" > ${log} 2>&1 &

    ((current_jobs += 1))

    if ((current_jobs >= max_jobs)); then
        wait -n
        ((current_jobs--))
    fi

done