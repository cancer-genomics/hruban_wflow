#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=6G
#$ -l h_vmem=6G
#$ -pe local 8
#$ -t 1-188
#$ -l cancergen

#-------
# Input
#--------------------------------------------------------------------------------------
manifest=../data/manifest.txt
outDir=../outDir
nProcesses=8
#--------------------------------------------------------------------------------------

module load conda_R/4.0.x

mkdir -p $outDir/plasma_obs

plasmaID=$(awk -F '\t' '{ if ($3=="plasma") print $2 }' $manifest | head -n $SGE_TASK_ID | tail -n 1)
plasmaAln=$(awk -F '\t' '{ if ($3=="plasma") print $5 }' $manifest | head -n $SGE_TASK_ID | tail -n 1)
patientID=$(awk -F '\t' '{ if ($3=="plasma") print $1 }' $manifest | head -n $SGE_TASK_ID | tail -n 1)
readLength=$(awk -F '\t' '{ if ($3=="plasma") print $4 }' $manifest | head -n $SGE_TASK_ID | tail -n 1)
tumorID=$(awk -F '\t' '{ if ($1 == "'$patientID'" && $3=="tumor") {print $2} }' $manifest)

strelkaSbsVcf=$outDir/strelka-somatic/$tumorID/results/variants/somatic.snvs.vcf.gz
outDir=$outDir/plasma_obs

mkdir -p $outDir/regions $outDir/pileup $outDir/obs $outDir/fragments $outDir/n_bases

Rscript ./i04-plasma-sbs-detection.R $strelkaSbsVcf $plasmaID $plasmaAln $outDir $nProcesses

