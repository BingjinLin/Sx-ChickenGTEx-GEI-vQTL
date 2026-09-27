#!/bin/bash
tissue=$1

dir_nominal="/ossfs/GxE/Etime/${tissue}"
dir_output="/data/Etime/${tissue}/Random_pairs"

mkdir -p ${dir_output}

### 1. output random nominal SNP-gene pairs information from nominal results to a .txt file
# (colnames: phenotype_id,variant_id,chr,pos)
# random select 100000 SNP-gene pairs in each tissue

  Etime=(Etime1 Etime2 Etime3)
  for time in "${Etime[@]}"
  do
  if [ -r "${dir_nominal}/${time}/${tissue}_${time}.cis_qtl.txt.gz" ]; then
      zcat ${dir_nominal}/${time}/${tissue}_${time}.cis_qtl_pairs.*.txt.gz | cut -f1,2 | grep -v '^pheno_id' | shuf -n 100000 > ${dir_output}/${tissue}_${time}.nominal_paires.shuf.txt
  fi
  done
   sort -u ${dir_output}/*.nominal_paires.shuf.txt > ${dir_output}/nominal_paires.shuf.txt


###Rscript===#
source /soft/miniconda3/etc/profile.d/conda.sh
conda activate myenv
Rscript /soft/zhh/vQTL/GXE/Etime/MashR/2-2Prepare-MashR-RandomPairs.R ${tissue}
bgzip "${dir_output}/nominal_pairs.combined_signifpairs.txt"
#> output file: nominal_pairs.combined_signifpairs.txt.gz



### 2. extract pairs from nominal results for each tissue#
source /soft/miniconda3/etc/profile.d/conda.sh
conda activate mashr

  Etime=(Etime1 Etime2 Etime3)
  for time in "${Etime[@]}"
  do
  if [ -r "${dir_nominal}/${time}/${tissue}_${time}.cis_qtl.txt.gz" ]; then
    # extract_pairs
    nominal_files=(`ls ${dir_nominal}/${time}/${tissue}_${time}.cis_qtl_pairs.*.txt.gz`)
       rm -f ${dir_output}/${tissue}_${time}.nominal_files.txt
       for l in ${nominal_files[*]}
        do
        {
            echo ${l} >> ${dir_output}/${tissue}_${time}.nominal_files.txt
        }
        done

    python3 /soft/zhh/vQTL/GXE/Etime/MashR/extract_pairs_lbj.py ${dir_output}/${tissue}_${time}.nominal_files.txt ${dir_output}/nominal_pairs.combined_signifpairs.txt.gz ${time}_nominal_pairs -o ${dir_output}
    #> output file: ${dir_output}/*.extracted_pairs.txt.gz
  fi
  done

### 3. prepare 1M random SNP-gene pairs for MashR
ls ${dir_output}/*.extracted_pairs.txt.gz > ${dir_output}/random_pairs_files.txt
# MashR format file (z-score)
python3 /soft/zhh/vQTL/GXE/Etime/MashR/mashr_prepare_input.py ${dir_output}/random_pairs_files.txt random_pairs -o ${dir_output} --only_zscore --dropna --subset 99999 --seed 9823
