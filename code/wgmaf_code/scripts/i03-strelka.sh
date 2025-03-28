#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=4G
#$ -l h_vmem=4G
#$ -pe local 8
#$ -t 1-34
#$ -l cancergen

module load python/2.7.9

#-------
# Input
#-----------------------------------------------------------------------------------------------
manifest=../data/manifest.txt
outDir=../outDir
refFasta= # Path to hg19 reference genome fasta file (not included). Contains chr1-22,X,Y,M.
callRegions=../data/call-regions-hg19.bed.gz # Included in ../data
STRELKA_INSTALL_PATH= # Path to Strelka installation (not included).
#-----------------------------------------------------------------------------------------------

mkdir -p $outDir/strelka-germline
mkdir -p $outDir/strelka-somatic

patientID=$(awk -F '\t' '{ if ($3=="tumor") print $1 }' $manifest | uniq | head -n $SGE_TASK_ID | tail -n 1)
tumorID=$(awk -F '\t' '{ if ($1 == "'$patientID'" && $3=="tumor") {print $2} }' $manifest)
tumorAln=$(awk -F '\t' '{ if ($1 == "'$patientID'" && $3=="tumor") {print $5} }' $manifest)
normalID=$(awk -F '\t' '{ if ($1 == "'$patientID'" && $3=="normal") {print $2} }' $manifest)
normalAln=$(awk -F '\t' '{ if ($1 == "'$patientID'" && $3=="normal") {print $5} }' $manifest)

#-------------------------
# Call somatic variants
#------------------------------------------------------------------------
python $STRELKA_INSTALL_PATH/bin/configureStrelkaSomaticWorkflow.py \
      --normalBam $normalAln \
      --tumorBam $tumorAln \
      --referenceFasta $refFasta \
      --runDir $outDir/strelka-somatic/$tumorID \
      --callRegions $callRegions

python $outDir/strelka-somatic/$tumorID/runWorkflow.py -m local -j 8
#------------------------------------------------------------------------


