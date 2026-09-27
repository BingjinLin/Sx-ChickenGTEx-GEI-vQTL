#!/bin/bash
date
dir_nominal="/juice_data/zhh_independent/ZHH_update/vQTL/res"
dir_output1="/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/top_pairs"

#mkdir -p ${dir_output1}

### 1. output all top SNP-gene pairs information from permutation results to a .txt file
echo " process step1 ......"
Rscript /soft/zhh/vQTL/tissue_sharing/1.1.combine_signif_pairs.R
echo " step1 job has finished......"

### 2. extract top pairs from nominal results for each tissue
echo " process step2 ......"
bash /soft/zhh/vQTL/tissue_sharing/1.2Prepare-MashR-StrongPairs.sh ${dir_nominal} ${dir_output1}
echo " step2 job has finished......"

### 3. prepare strong SNP-gene pairs for MashR
echo " process step3 ......"
bash /soft/zhh/vQTL/tissue_sharing/1.3Prepare-MashR-StrongPairs.sh ${dir_nominal} ${dir_output1}
echo " step3 job has finished......"
