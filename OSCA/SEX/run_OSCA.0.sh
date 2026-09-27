#!/bin/bash
set -e

TOTAL_TASKS=32
MAX_PARALLEL=8

#mkdir -p logs

tissue_list=($(cat /soft/zhh/vQTL/OSCA/tissue.list))
len=${#tissue_list[@]}
tissue_list=("${tissue_list[@]:0:$((len-2))}")

for tis in "${tissue_list[@]}"
do
  SEX=(male female)

  for Sx in "${SEX[@]}"
  do
      echo "Launching task for $tis $Sx..."

      bash /soft/zhh/vQTL/OSCA/SEX/OSCA_pre_1.sh "$tis" "$Sx" \
        > "./logs/run.$tis.$Sx.log" 2>&1 &

      if (( $(jobs -r | wc -l) >= MAX_PARALLEL ))
      then
        wait -n
      fi
  done
done

wait
echo "All tasks completed."