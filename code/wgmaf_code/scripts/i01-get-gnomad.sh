#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=100G
#$ -l h_vmem=100G
#$ -l cancergen

module load conda_R/4.0.x

#--------
# Input
#-------------------
tmpDir=../data/temp
outDir=../data
#-------------------

mkdir -p $tmpDir $outDir

wget http://hgdownload-euro.soe.ucsc.edu/gbdb/hg38/gnomAD/vcf/gnomad.genomes.r3.0.sites.vcf.gz -O $tmpDir/gnomad.genomes.r3.0.sites.vcf.gz
wget http://hgdownload-euro.soe.ucsc.edu/gbdb/hg38/gnomAD/vcf/gnomad.genomes.r3.0.sites.vcf.gz.tbi -O $tmpDir/gnomad.genomes.r3.0.sites.vcf.gz.tbi

Rscript ./i01-get-gnomad.R $tmpDir $outDir
