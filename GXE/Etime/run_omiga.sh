#!/bin/bash
set -e

MAX_PARALLEL=3

mkdir -p logs

tissue_list=($(cat /soft/zhh/vQTL/OSCA/tissue.list))
tis="Testis"
#for tis in "${tissue_list[@]}"
#do
  Etime=(1 2 3)

  for time in "${Etime[@]}"
  do

      bash Omiga_inte.1.sh "$tis" "$time" \
        > "./logs/run.$tis.$time.log" 2>&1 &

      if (( $(jobs -r | wc -l) >= MAX_PARALLEL ))
      then
        wait -n
      fi
  done
#done

wait
echo "All tasks completed."