#!/bin/bash
date
tissue=$1

dir_nominal="/ossfs/GxE/Etime/${tissue}"
dir_output="/data/Etime/${tissue}/strong_pairs"

mkdir -p ${dir_output}

### 1. output all top SNP-gene pairs information from permutation results to a .txt file
#echo "Process step1.1 analysis on tissue: ${tissue}"
#bash /soft/zhh/vQTL/GXE/Etime/MashR/1-1Prepare-MashR-StrongPairs.sh ${dir_nominal} ${dir_output}

### 2. extract top pairs from nominal results for each tissue
echo "Process step1.2 analysis on tissue: ${tissue}"
bash /soft/zhh/vQTL/GXE/Etime/MashR/1-2Prepare-MashR-StrongPairs.sh ${dir_nominal} ${dir_output} ${tissue}

### 3. prepare strong SNP-gene pairs for MashR
echo "Process step1.3 analysis on tissue: ${tissue}"
bash /soft/zhh/vQTL/GXE/Etime/MashR/1-3Prepare-MashR-StrongPairs.sh ${dir_nominal} ${dir_output}

date
