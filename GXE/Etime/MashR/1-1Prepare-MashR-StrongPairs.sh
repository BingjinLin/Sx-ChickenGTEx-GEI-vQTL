#!/bin/bash
dir_nominal=${1}
dir_output=${2}


############################################################################################
#   1: combine_signif_pairs_lbj.py   2: extract_pairs_lbj.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 1. output all top SNP-gene pairs information from permutation results to a .txt file
# (colnames: phenotype_id,variant_id,chr,pos)
# file list of permutation results 
rm -f ${dir_output}/permutation_files.txt
# permutation results
perm_files=($(find "${dir_nominal}" -type f -name "*.cis_qtl.txt.gz")) # permutation results
for l in ${perm_files[*]}
do
{
    echo ${l} >> ${dir_output}/permutation_files.txt
}
done

python3 /soft/zhh/vQTL/GXE/Etime/MashR/combine_signif_pairs_lbj.py ${dir_output}/permutation_files.txt strong_pairs -o ${dir_output}


