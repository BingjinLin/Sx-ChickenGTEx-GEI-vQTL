#!/bin/bash

TOTAL_TASKS=39
MAX_PARALLEL=16

#mkdir -p logs

tissue_list=($(cat /soft/zhh/vQTL/OSCA/tissue.list))
len=${#tissue_list[@]}
tissue_list=("${tissue_list[@]:0:$((len-2))}")

for tis in "${tissue_list[@]}"
do
  SEX=(male female)

  for Sx in "${SEX[@]}"
  do

    for ((i=1; i<=TOTAL_TASKS; i++))
    do

      while [ "$(jobs -rp | wc -l)" -ge "$MAX_PARALLEL" ]
      do
          wait -n
      done

      echo "Launching task for $tis $Sx , chr=$i..."

(
    logfile="./logs/run.${tis}.${Sx}.${i}.log"

    if bash OSCA_vQTL_chr.2.sh "$tis" "$Sx" "$i" \
        > "$logfile" 2>&1
    then
        rm -f "$logfile"
    else
        echo "Task failed: $tis $Sx chr=$i"
    fi
) &
    done
  done
done

wait
echo "All tasks completed."