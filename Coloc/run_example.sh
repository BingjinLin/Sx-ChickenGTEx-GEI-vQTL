#!/bin/bash
set -e

MAX_PARALLEL=3

source /soft/miniconda3/etc/profile.d/conda.sh
conda activate myenv
#script="/soft/zhh/vQTL/Coloc/example_plot_tGWAS_var.R"
script="/soft/zhh/vQTL/Coloc/example.plot_ieQTL_vQTL_mGWAS.R"

tail -n +2 /juice_data/zhh_independent/ZHH_update/vQTL/res/Plot/colocation/example/vQTL_mGWAS_ieQTL/task.list | \
while read tissue gene chr trait loci type interaction
do
    echo $tissue $gene $chr $trait $loci

log="./logs/${tissue}_${gene}_${loci}.plot.log"
nohup Rscript ${script} "${tissue}" "${gene}" "${chr}" "${trait}"  "${loci}" > ${log} 2>&1 || echo "Error in ${tissue}_${gene}_${loci}" &

nohup Rscript ${script} "${tissue}" "${gene}" "${chr}" "${trait}"  "${loci}" "${type}" "${interaction}" > ${log} 2>&1 || echo "Error in ${tissue}_${gene}_${loci}" &

      if (( $(jobs -r | wc -l) >= MAX_PARALLEL ))
      then
        wait -n
      fi

done
wait
echo "All tasks completed."