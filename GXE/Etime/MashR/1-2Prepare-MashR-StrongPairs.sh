#!/bin/bash

dir_nominal=${1}
dir_output=${2}
tissue=${3}
############################################################################################
#   1: combine_signif_pairs_lbj.py   2: extract_pairs_lbj.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 2. extract top pairs from nominal results for each tissue
Etimes=(Etime1 Etime2 Etime3)

# echo $nFile # 0..32
for name in "${Etimes[@]}"
do
{
  if [ -r "${dir_nominal}/${name}/${tissue}_${name}.cis_qtl.txt.gz" ]; then
    nominal_files=(`ls ${dir_nominal}/${name}/${tissue}_${name}.cis_qtl_pairs.*.txt.gz`)
    rm -f ${dir_output}/${name}.nominal_files.txt
    for l in ${nominal_files[*]}
    do
    {
        echo ${l} >> ${dir_output}/${name}.nominal_files.txt
    }
    done
    # extract_pairs
    python3 /soft/zhh/vQTL/GXE/Etime/MashR/extract_pairs_lbj.py ${dir_output}/${name}.nominal_files.txt ${dir_output}/strong_pairs.combined_signifpairs.txt.gz ${name} -o ${dir_output}
    ### output file: *.extracted_pairs.txt.gz
    rm -f ${dir_output}/${name}.nominal_files.txt
  fi
} &
done
wait

