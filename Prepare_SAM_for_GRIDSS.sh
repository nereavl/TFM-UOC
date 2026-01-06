#!/bin/bash

# Processing SAM/BAM files from the mapping (alignment) of OF74 reads using samtools
STRAIN=$1
BAM_ALINGMENT=$2
HOLOGENOME=$3
DSUB_CHR=$4  # RefSeq chromosome code (e.g., NC_048534.1)

mkdir -p seqs/3_SAM_for_GRIDSS

# copy and index chosen BAM alignment
echo $"Copying and indexing BAM alignment"
cp $BAM_ALINGMENT ./seqs/3_SAM_for_GRIDSS/${STRAIN}_mapped.bam
samtools index ./seqs/3_SAM_for_GRIDSS/${STRAIN}_mapped.bam

# copy and index chosen BAM hologenome file
echo $"Copying and indexing hologenome"
cp $HOLOGENOME ./seqs/3_SAM_for_GRIDSS/hologenome.fa
samtools faidx ./seqs/3_SAM_for_GRIDSS/hologenome.fa --write-index

for ii in $DSUB_CHR; do
    echo $'\n\n'"Detecting breakpoints in chromosome: $ii"
	awk '{print $1, 0, $2}' OFS="\t" ./seqs/3_SAM_for_GRIDSS/hologenome.fa.fai > ./seqs/3_SAM_for_GRIDSS/full.bed
	grep -v "Drosophila_subobscura_$ii" ./seqs/3_SAM_for_GRIDSS/full.bed > ./seqs/3_SAM_for_GRIDSS/$ii.bed
done
