#!/bin/bash
set -e

MAX_PARALLEL=3

mkdir -p logs

script="/soft/zhh/vQTL/GXE/Etime/MashR/01.prepare_MashR_Strong.sh"
source /soft/miniconda3/etc/profile.d/conda.sh
conda activate mashr

#tissue_list=($(cat /soft/zhh/vQTL/OSCA/tissue.list))
#tissue_list=("${tissue_list[@]:1:32}")
tissue_list=(Testis Ovary)
for tis in "${tissue_list[@]}"
do
  echo "Process step1 analysis on tissue: ${tis}"
  bash ${script} "$tis" \
        > "./logs/${tis}_pre1.log" 2>&1 &

      if (( $(jobs -r | wc -l) >= MAX_PARALLEL ))
      then
        wait -n
      fi

done

wait
echo "All tasks completed."