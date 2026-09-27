#!/bin/bash
set -eu
# 控制并发任务数
max_jobs=4
current_jobs=0

source /soft/miniconda3/etc/profile.d/conda.sh
conda activate TWAS

 for tis in `cat /ossfs/OmiGA/tiss.txt`
 do

Otfile="/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/Random_pairs/${tis}_nominal_pairs.extracted_pairs.txt.gz"
if [ ! -f "$Otfile" ]; then

    echo "----------------------------------------"
    echo "Submitting job for: tissue=$tis"

    log_file="/soft/zhh/vQTL/tissue_sharing/logs/${tis}_sharing.log"
    Rscript /soft/zhh/vQTL/tissue_sharing/extract_pairs_lbj.R "${tis}" > "${log_file}" 2>&1 &

    ((current_jobs += 1))

    if ((current_jobs >= max_jobs)); then
        wait -n  # 等待任意一个后台任务完成
        ((current_jobs--))  # 减少一个任务计数
    fi
fi
done

  echo "All jobs have been submitted. Waiting for them to finish..."
  wait # 等待所有后台进程结束
  echo "All jobs have been finished."
  
