#!/bin/bash

source /soft/miniconda3/etc/profile.d/conda.sh
conda activate myenv

dir_input="/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing"
strong_file="${dir_input}/top_pairs/top_pairs.MashR_input.txt.gz"

random_file="${dir_input}/Random_pairs/random_pairs.MashR_input.txt.gz"
dir_output="${dir_input}/result"
mkdir -p ${dir_output}

Rscript /soft/zhh/vQTL/tissue_sharing/run_MashR.R ${strong_file} ${random_file} 0 ${dir_output}
