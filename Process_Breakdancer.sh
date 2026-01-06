#!/bin/bash

############################################
# Script: run_breakdancer.sh
# Description: Detects chromosomal breakpoints
#              using BreakDancer in
# Author: Mercé Merayo & Kenia Delgado
# Adaptation: Nerea Vallejo
############################################

# ----------------------------
# Positional arguments
# ----------------------------
STRAIN=$1
BAM_ALIGNMENT=$2
DSUB_CHR=$3

# ----------------------------
# Prepare folders
# ----------------------------
echo $"Preparing folders"
mkdir -p seqs/3_SAM_for_breakdancer
mkdir -p ./seqs/4_BreakDancer
mkdir -p ./seqs/5_Clean_BreakDancer

# ----------------------------
# Pre process files
# ----------------------------
echo $"Copying and indexing BAM alignment"
cp $BAM_ALIGNMENT "./seqs/3_SAM_for_breakdancer/${STRAIN}_mapped.bam"
samtools index ./seqs/3_SAM_for_breakdancer/${STRAIN}_mapped.bam

# ----------------------------
# Run & Clean Breakdancer
# ----------------------------
echo "Generating configuration file for Breakdancer..."
bam2cfg.pl -g -h ./seqs/3_SAM_for_breakdancer/${STRAIN}_mapped.bam > "./seqs/4_BreakDancer/${STRAIN}_BP_bwamem2.cfg"

for ii in $DSUB_CHR; do
    echo $'\n\n'"Detecting breakpoints in chromosome: $ii"
	breakdancer-max -o "Drosophila_subobscura_$ii" "./seqs/4_BreakDancer/${STRAIN}_BP_bwamem2.cfg" -y 0 -q -1 -g "./seqs/4_BreakDancer/${STRAIN}_breakdancer_$ii.bed" > "./seqs/4_BreakDancer/${STRAIN}_breakdancer_$ii.txt"
	echo $'\n\n'"Cleaning outputs of chromosome: $ii"
	Rscript ./scripts_to_publish/Clean_Breakdancer.R $STRAIN $ii 1000000 0 > ./seqs/5_Clean_BreakDancer/${ii}_log.txt
done
