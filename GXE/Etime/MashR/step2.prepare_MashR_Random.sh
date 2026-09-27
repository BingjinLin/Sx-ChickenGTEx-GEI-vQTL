#!/bin/bash
date
set -e 

chr=$1

dir_output="/juice_data/splicing_temp/sQTLmapping/tissue_sharing/Random_pairs"
dir_nominal="/juice_data/splicing_temp/sQTLmapping"


mkdir -p ${dir_output}

### 1. output all top SNP-gene pairs information from permutation results to a .txt file
###768GB 内存服务器,chunk 不宜太小， 1024个 或者 所有数据放入1个chunk 即可
bash /soft/zhh/MashR/2-1Prepare-MashR-RandomPairs.sh ${dir_nominal} ${dir_output} $chr
#upload $dir_output $ossfs/tissue_sharing/MashR/result2/

### 2. extract top pairs from nominal results for each tissue
#256GB 内存运行一个组织 需要大内存服务器。
#bash /ossfs/tissue_sharing/MashR/2-2Prepare-MashR-RandomPairs_ycg.sh ${dir_nominal} ${dir_output}
#upload $dir_output $/ossfs/tissue_sharing/MashR/result/Random_pairs/

### 3. prepare strong SNP-gene pairs for MashR
#bash /ossfs/tissue_sharing/MashR/2-3Prepare-MashR-RandomPairs.sh ${dir_output}

date
