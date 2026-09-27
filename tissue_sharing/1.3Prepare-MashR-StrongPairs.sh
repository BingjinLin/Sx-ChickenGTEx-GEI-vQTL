#!/bin/bash

dir_nominal=${1}
dir_output=${2}

############################################################################################
#   1: combine_signif_pairs_lbj.py   2: extract_pairs_lbj.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 3. prepare four SNP-gene pairs for MashR
top_pairs_files=(`ls ${dir_output}/*_strong_pairs.extracted_pairs.txt.gz`)
rm -f ${dir_output}/top_pairs_files.txt
for l in ${top_pairs_files[*]}
do
{
    echo ${l} >> ${dir_output}/top_pairs_files.txt
}
done
# MashR format file (z-score)
python3 /soft/zhh/vQTL/tissue_sharing/mashr_prepare_input.py ${dir_output}/top_pairs_files.txt top_pairs -o ${dir_output} --only_zscore
zcat ${dir_output}/top_pairs.MashR_input.txt.gz | sed -e 's/_zval//g' | gzip > ${dir_output}/top_pairs.temp.txt.gz
rm -f ${dir_output}/top_pairs.MashR_input.txt.gz
mv ${dir_output}/top_pairs.temp.txt.gz ${dir_output}/top_pairs.MashR_input.txt.gz



