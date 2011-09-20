#!/bin/sh -w
#$ -S /bin/sh
#$ -cwd
#$ -l mem=8G,time=24::

##Script to run vcfCodingSnps.v1.5 which annotates files in vcf format. Other required inputs are REference Genome and Gene list
## Inputs:
##### 1. Input VCF file (full path)
##### 2. Name of output directory
##### 3. Name of run i.e. prefix for .out .log files


##vcfCodingSnps installation directory
VCF_DIR="/ifs/scratch/c2b2/ip_lab/sz2317/privacy/workingdir/vcfCodingSnps.v1.5"

## vcfCodingSnps program
VCFPROG=$VCF_DIR"/vcfCodingSnps.v1.5"

## REference genome and geneList
REF_GENE=$VCF_DIR"/referenceGenomes/genome.V36.fa"          #default used is B36 for human
GENE_LIST=$VCF_DIR"/geneLists/UCSCknownGene.B36.txt"		#default used in UCSC.B36 for human
## input file
PREFIX_INPUT=$1

## path to output directory
PREFIX_RUN=$VCF_DIR"/runs/"$2

## name of run
NAME=$3

echo "$VCFPROG -s $PREFIX_INPUT --refgenome $REF_GENE --genefile $GENE_LIST -o $PREFIX_RUN"/"$NAME".out" -l $PREFIX_RUN"/"$NAME".log" " >  $PREFIX_RUN"/log"

$VCFPROG -s $PREFIX_INPUT --refgenome $REF_GENE --genefile $GENE_LIST -o $PREFIX_RUN"/"$NAME".out" -l $PREFIX_RUN"/"$NAME".log"




~

