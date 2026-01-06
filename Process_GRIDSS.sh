#!/bin/bash

# Detect Structural Variants (SVs) and breakpoints using LUMPY
# The generated txt file is stored in a separate folder to improve result visibility.

# ----------------------------
# Positional arguments
# ----------------------------
STRAIN=$1
BAM_ALINGMENT=$2
HOLOGENOME=$3
DSUB_CHR=$4

# ----------------------------
# Prepare folders
# ----------------------------
echo $"Preparing folders"
mkdir -p ./seqs/3_SAM_for_GRIDSS
mkdir -p ./seqs/4_GRIDSS
mkdir -p ./seqs/5_Clean_GRIDSS


# copy and index chosen BAM alignment
echo $"Copying and indexing BAM alignment"
# cp $BAM_ALINGMENT ./seqs/3_SAM_for_GRIDSS/${STRAIN}_mapped.bam
# samtools index ./seqs/3_SAM_for_GRIDSS/${STRAIN}_mapped.bam

# copy and index chosen BAM hologenome file
echo $"Copying and indexing hologenome"
# cp $HOLOGENOME ./seqs/3_SAM_for_GRIDSS/hologenome.fa
# samtools faidx ./seqs/3_SAM_for_GRIDSS/hologenome.fa --write-index

for ii in $DSUB_CHR; do
   
    echo $'\n\n'"Detecting breakpoints in chromosome: $ii"
	awk '{print $1, 0, $2}' OFS="\t" ./seqs/3_SAM_for_GRIDSS/hologenome.fa.fai > ./seqs/3_SAM_for_GRIDSS/full.bed
	grep -v "Drosophila_subobscura_$ii" ./seqs/3_SAM_for_GRIDSS/full.bed > ./seqs/3_SAM_for_GRIDSS/$ii.bed

    echo $'\n\n'"Detecting breakpoints in chromosome: $ii"
	mkdir -p ./seqs/4_GRIDSS/$ii
	awk '{print $1, 0, $2}' OFS="\t" ./seqs/3_SAM_for_GRIDSS/hologenome.fa.fai > ./seqs/3_SAM_for_GRIDSS/full.bed
	grep -v "Drosophila_subobscura_$ii" ./seqs/3_SAM_for_GRIDSS/full.bed > ./seqs/3_SAM_for_GRIDSS/$ii.bed
	
	echo $"Running GRIDSS"
	gridss -r ./seqs/3_SAM_for_GRIDSS/hologenome.fa -b ./seqs/3_SAM_for_GRIDSS/$ii.bed ./seqs/3_SAM_for_GRIDSS/${STRAIN}_mapped.bam -j ~/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar -o ./seqs/4_GRIDSS/$ii/$ii.vcf -a ./seqs/4_GRIDSS/$ii/$ii.bam -w ./seqs/4_GRIDSS
	
	samtools index ./seqs/4_GRIDSS/$ii/$ii.bam
	
	echo $'\n\n'"Processing SV calling in chromosome $ii"
    Rscript ./scripts_to_publish/Clean_GRIDSS.R $STRAIN $ii 1000000 0 > ./seqs/5_Clean_GRIDSS/${ii}_log.txt
	
done