#!/bin/bash
set -e
set -x

task_file="/juice_data/zhh_independent/ZHH_update/vQTL/res/vqtl_GWAS_coloc.plot_task.txt"
max_jobs=2
current_jobs=0

tail -n +2 ${task_file} | while read tissue gene chrom trait bin
do
    nohup Rscript example_plot_pre.R $tissue $gene $chrom $trait $bin \
    > ./plot_logs/run_${tissue}_${gene}.log 2>&1 &

    ((current_jobs += 1))

    if ((current_jobs >= max_jobs)); then
        wait -n
        ((current_jobs--))
    fi

done