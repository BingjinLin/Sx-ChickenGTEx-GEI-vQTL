#!/bin/bash
set -e
tissue=$1 
trait=$2

omiga="/soft/OmiGA-build-v1.1.3-beta.2+250831/bin/omiga"
Datadir="/juice_data/zhh_independent/ZHH_update"
covar="${Datadir}/covariates/${tissue}.auto_covar"
tmp_dir="${Datadir}/GxE/veQTL/tmp"
expr="${tmp_dir}/${tissue}_${trait}.expr_tmm_inv.bed"
geno="${tmp_dir}/${tissue}_${trait}.geno"
interactions_file="${tmp_dir}/${tissue}_${trait}.inte.txt"
outputdir="/juice_data/zhh_independent/ZHH_update/GxE/veQTL/${tissue}/${trait}"
mkdir -p $outputdir
prefix=${tissue}_${trait}

/usr/bin/time ${omiga} --mode cis_interaction --qtl-map-model a \
           --genotype ${geno} --phenotype ${expr} --covariates ${covar} \
           --normalized-interaction --interaction ${interactions_file} --interaction-name ${trait} \
           --prefix ${prefix} --output-dir ${outputdir}
rm -f ${expr}
rm -f ${geno}.bed
rm -f ${geno}.bim
rm -f ${geno}.fam
rm -f ${geno}.nosex
rm -f ${interactions_file}

echo "Omiga process has been done!"