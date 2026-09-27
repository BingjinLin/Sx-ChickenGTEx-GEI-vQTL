#!/bin/bash
set -e

MAX_PARALLEL=3

source /soft/miniconda3/etc/profile.d/conda.sh
conda activate myenv
script="/soft/zhh/vQTL/Coloc/coloc_evQTL_Etime.R"

tissue_list=($(cat /soft/zhh/vQTL/OSCA/tissue.list))
tissue_list=("${tissue_list[@]:29:32}")
for tis in "${tissue_list[@]}"
do
  Etime=(Etime1 Etime2 Etime3)

  for time in "${Etime[@]}"
  do

echo " processing coloc analysis on ${tis} : ${time}"
log="./logs/${tis}_${time}.coloc.log"
nohup Rscript ${script} "${tis}" "${time}" > ${log} 2>&1 &

      if (( $(jobs -r | wc -l) >= MAX_PARALLEL ))
      then
        wait -n
      fi
  done
done

wait
echo "All tasks completed."