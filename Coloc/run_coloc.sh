#!/bin/bash
set -e

tissue_list=($(cat /soft/zhh/vQTL/OSCA/tissue.list))
#tissue_list=("${tissue_list[@]:29:32}")
#tissue_list=("Whole.blood" "Testis")
#tissue_list=("Crop" "Duodenum" "Fat" "Gizzard" "Heart" "Hypothalamus" "Ileum" "Jejunum" "Kidney" "Leg.muscle" "Whole.blood" "Lung") #sex
#tissue_list=("Adrenals" "Duodenum" "Gizzard" "Heart" "Hypothalamus" "Ileum" "Jejunum" "Kidney" "Leg.muscle" "Whole.blood" "Lung" "Ovary" "Pancreas" "Retina" "Skin" "Testis" "Trachea") #CELL

max_jobs=4
current_jobs=0

for tis in "${tissue_list[@]}"
#cat /soft/zhh/vQTL/OSCA/tiss_cell.list | while read tis cell

do

echo "Launching task for $tis..."
  Rscript /soft/zhh/vQTL/Coloc/coloc_evQTL_1.R "$tis" > "./logs/run.$tis.ecoloc.log" 2>&1 &
  #Rscript /soft/zhh/vQTL/OSCA/coloc_ievQTL.R "$tis" > "./logs/run.$tis.iecoloc.log" 2>&1 &
  #Rscript /soft/zhh/vQTL/OSCA/coloc_vQTL_sex.R "$tis" > "./logs/run.$tis.scoloc.log" 2>&1 &

    ((current_jobs += 1))

    if ((current_jobs >= max_jobs)); then
        wait -n
         ((current_jobs -= 1))
    fi

done

wait
echo "All tasks completed."
