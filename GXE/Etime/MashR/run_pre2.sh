#!/bin/bash
set -e

MAX_PARALLEL=3

mkdir -p logs

script="/soft/zhh/vQTL/GXE/Etime/MashR/02.Prepare-MashR-RandomPairs.sh"

#tissue_list=($(cat /soft/zhh/vQTL/OSCA/tissue.list))
#tissue_list=("${tissue_list[@]:3:32}")
tissue_list=(Testis Ovary)

for tis in "${tissue_list[@]}"
do
  echo "Process step2 analysis on tissue: ${tis}"
  bash ${script} "$tis" \
        > "./logs/${tis}_pre2.log" 2>&1 &

      if (( $(jobs -r | wc -l) >= MAX_PARALLEL ))
      then
        wait -n
      fi

done

wait
echo "All tasks completed."