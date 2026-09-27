#!/bin/bash
set -e

MAX_PARALLEL=3

source /soft/miniconda3/etc/profile.d/conda.sh
conda activate myenv
script="/soft/zhh/vQTL/Coloc/example_plot_tGWAS.R"
Task_file="/juice_data/zhh_independent/ZHH_update/vQTL/res/Plot/high_vcoloc.tOSCA.task.list"

tail -n +2 "${Task_file}" | while read tissue gene chr trait loci
do
    echo $tissue $gene $chr $trait $loci

log="./logs/example/${tissue}_${gene}_${chr}_${loci}.plot.log"
(nohup Rscript ${script} "${tissue}" "${gene}" "${chr}" "${trait}" "${loci}" > ${log} 2>&1 || echo "Error in ${tissue}_${gene}_${chr}_${loci}") &
      if (( $(jobs -r | wc -l) >= ${MAX_PARALLEL} ))
      then
        wait -n
      fi

done
wait
echo "All tasks completed."