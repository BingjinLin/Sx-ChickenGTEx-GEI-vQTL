#!/bin/bash
set -e

MAX_PARALLEL=3

source /soft/miniconda3/etc/profile.d/conda.sh
conda activate myenv
script="/soft/zhh/vQTL/Coloc/coloc.ieQTL_G.R"

tissue_list=($(cat /soft/zhh/vQTL/OSCA/tissue.list))
for tis in "${tissue_list[@]}"
do

echo " processing coloc analysis on ${tis} :"
log="./logs/${tis}_Growth.coloc.log"
nohup Rscript ${script} "${tis}" > ${log} 2>&1 &

      if (( $(jobs -r | wc -l) >= MAX_PARALLEL ))
      then
        wait -n
      fi

done

wait
echo "All tasks completed."