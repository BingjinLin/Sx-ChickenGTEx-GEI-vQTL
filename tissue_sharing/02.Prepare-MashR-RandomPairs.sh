#!/bin/bash

dir_output="/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/Random_pairs"
dir_nominal="/juice_data/zhh_independent/ZHH_update/vQTL/res"

mkdir -p ${dir_output}

### 1. output random nominal SNP-gene pairs information from nominal results to a .txt file
# (colnames: phenotype_id,variant_id,chr,pos)
# random select 100000 SNP-gene pairs in each tissue

#==for tis in `cat /ossfs/OmiGA/tiss.txt`
#==do
#==       zcat ${dir_nominal}/${tis}/zscore/${tis}.vqtl_chr*.all.query.gz | cut -f10,1 | grep -v '^Gene' | shuf -n 100000 > ${dir_output}/${tis}.nominal_paires.shuf.txt
#==done

#==sort -u ${dir_output}/*.nominal_paires.shuf.txt > ${dir_output}/nominal_paires.shuf.txt
#===Rscript===#
#library("data.table")
#==data<-fread("/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/Random_pairs/nominal_paires.shuf.txt",head=F)
#==info<-fread("/ossfs/tissue_sharing/MashR/GRCg7b_SNP.txt")
#==df<-data.frame(Chr=info$V1[match(data$V1,info$V2)],SNP=data$V1, BP=info$V3[match(data$V1,info$V2)], Gene=data$V2)
#==write.table(df, "/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/Random_pairs/nominal_pairs.combined_signifpairs.txt",sep="\t",quote=F,row.names=F)
bgzip ${dir_output}/nominal_pairs.combined_signifpairs.txt
#> output file: nominal_pairs.combined_signifpairs.txt.gz

### 2. extract pairs from nominal results for each tissue#
#==bash /soft/zhh/vQTL/tissue_sharing/2.2Prepare-MashR-RandomPairs.sh

### 3. prepare 1M random SNP-gene pairs for MashR
ls ${dir_output}/*_nominal_pairs.extracted_pairs.txt.gz > ${dir_output}/random_pairs_files.txt
# MashR format file (z-score)
python3 /soft/zhh/vQTL/tissue_sharing/mashr_prepare_input.py ${dir_output}/random_pairs_files.txt random_pairs -o ${dir_output} --only_zscore --dropna --subset 1000000 --seed 9823
