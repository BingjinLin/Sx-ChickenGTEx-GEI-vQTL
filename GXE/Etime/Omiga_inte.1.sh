#!/bin/bash
set -e
tissue=$1 
Etime=$2

omiga="/soft/OmiGA-build-v1.1.3-beta.2+250831/bin/omiga"

covar="/data/bed/${tissue}_Etime${Etime}_known.covar.list"
expr="/data/bed/${tissue}_Etime${Etime}.expr_tmm_inv.bed.gz"
geno="/juice_data/zhh_independent/ZHH_update/Genotype/${tissue}/${tissue}_Etime${Etime}.autosomal.geno"

outputdir="/data/GxE/${tissue}/Etime${Etime}"
mkdir -p $outputdir
prefix=${tissue}_Etime${Etime}

if [ -f "${expr}" ] && [ -f "${covar}" ]; then

    echo "Starting analysis with covariates..."

/usr/bin/time ${omiga} --mode cis --qtl-map-model a+A \
           --genotype ${geno} --phenotype ${expr}  \
           --geno-pc-covar 3 --dprop-pc-covar 0.001 --rm-collinear-covar 0.95 \
           --covariates ${covar} --dcovar-name sex \
           --maf-threshold 0.05 \
           --output-dir ${outputdir} --prefix ${prefix}
else
if [ -f "${expr}" ]; then
    echo "Starting analysis without covariates (covar file missing)..."

/usr/bin/time ${omiga} --mode cis --qtl-map-model a+A \
           --genotype ${geno} --phenotype ${expr}  \
           --geno-pc-covar 3 --dprop-pc-covar 0.001 --rm-collinear-covar 0.95 \
           --maf-threshold 0.05 \
           --output-dir ${outputdir} --prefix ${prefix}
else
    echo "There is no combination analysis..."
fi
fi

echo "Omiga process has been done!"