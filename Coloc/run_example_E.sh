#!/bin/bash
set -e

MAX_PARALLEL=3

source /soft/miniconda3/etc/profile.d/conda.sh
conda activate myenv
script="/soft/zhh/vQTL/Coloc/example_plot_Eievcoloc_high.R"

tail -n +3 Etime_example.list | while read tissue gene chr time
do
    echo $tissue $gene $chr

log="./logs/${tissue}_${gene}_${chr}_${time}.plot.log"
nohup Rscript ${script} "${tissue}" "${gene}" "${chr}" "${time}" > ${log} 2>&1 || echo "Error in ${tissue}_${gene}_${chr}" &
      if (( $(jobs -r | wc -l) >= MAX_PARALLEL ))
      then
        wait -n
      fi

done
wait
echo "All tasks completed."