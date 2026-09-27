#!/bin/bash
set -e
set -x
#trait=$1
tail -n +2 /soft/zhh/GWAS/colocation/trait_info.txt | while IFS=$'\t' read -r trait _; do

echo " processing coloc analysis on $trait"
for i in $(seq 1 4 32); do
    end=$((i+3))

    echo "bash /soft/zhh/vQTL/Coloc/run_coloc_1.sh $i $end $trait" >> /soft/zhh/vQTL/Coloc/New_bash/All_coloc.vQTL.tGWAS.sh
done

done
