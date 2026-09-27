#!/bin/bash
set -e

MAX_PARALLEL=3

source /soft/miniconda3/etc/profile.d/conda.sh
conda activate myenv

dir_input="/ossfs/GxE/Etime/MashR"

tissue_list=($(cat /soft/zhh/vQTL/OSCA/tissue.list))
tissue_list=("${tissue_list[@]:1:32}")

for tissue in "${tissue_list[@]}"
do
  echo "Process MashR analysis on tissue: ${tissue}"

strong_file="${dir_input}/${tissue}/strong_pairs/strong_pairs.MashR_input.txt.gz"
random_file="${dir_input}/${tissue}/Random_pairs/random_pairs.MashR_input.txt.gz"
dir_output="/data/Etime/result/${tissue}"
mkdir -p ${dir_output}

Rscript /soft/zhh/vQTL/GXE/Etime/MashR/run_MashR.R ${strong_file} ${random_file} 0 ${dir_output} > "./logs/${tissue}_mashr.log" 2>&1 &


      if (( $(jobs -r | wc -l) >= MAX_PARALLEL ))
      then
        wait -n
      fi

done

wait
echo "All tasks completed."
