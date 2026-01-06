#!/bin/bash

# Detect Structural Variants (SVs) and breakpoints using LUMPY
# The generated txt file is stored in a separate folder to improve result visibility.

# ----------------------------
# Positional arguments
# ----------------------------
STRAIN=$1
BAM_ALINGMENT=$2
HOLOGENOME=$3
DSUB_CHR=$4  # RefSeq chromosome code (e.g., NC_048534.1)

# ----------------------------
# Prepare folders
# ----------------------------
echo $"Preparing folders"
mkdir -p ./seqs/3_SAM_for_Manta
mkdir -p ./seqs/4_Manta
mkdir -p ./seqs/5_Clean_Manta


# copy and index chosen BAM alignment
echo $"Copying and indexing BAM alignment"
cp $BAM_ALINGMENT ./seqs/3_SAM_for_Manta/${STRAIN}_mapped.bam
samtools index ./seqs/3_SAM_for_Manta/${STRAIN}_mapped.bam

# copy and index chosen BAM hologenome file
echo $"Copying and indexing hologenome"
cp $HOLOGENOME ./seqs/3_SAM_for_Manta/hologenome.fa
samtools faidx ./seqs/3_SAM_for_Manta/hologenome.fa --write-index


for ii in $DSUB_CHR; do
    echo $'\n\n'"PREPARING REGION BED FILES FOR MANTA, CHROMOSOME: $ii"
	
    awk -v chr="Drosophila_subobscura_$ii" '$1 == chr {print $1, 0, $2}' OFS="\t" ./seqs/3_SAM_for_Manta/hologenome.fa.fai > ./seqs/3_SAM_for_Manta/$ii.bed
	bgzip -f ./seqs/3_SAM_for_Manta/$ii.bed
    tabix -p bed ./seqs/3_SAM_for_Manta/$ii.bed.gz

    echo $'\n\n'"Detecting breakpoints in chromosome: $ii"
	
	echo $"Running Manta"
	mkdir -p seqs/4_Manta/$ii
	python2 /usr/bin/manta/bin/configManta.py --bam ./seqs/3_SAM_for_Manta/${STRAIN}_mapped.bam --referenceFasta ./seqs/3_SAM_for_Manta/hologenome.fa --runDir ./seqs/4_Manta/$ii --callRegions ./seqs/3_SAM_for_Manta/$ii.bed.gz
	
	./seqs/4_Manta/$ii/runWorkflow.py
	
	gzip -d -r seqs/4_Manta/$ii/results/variants/diploidSV.vcf
	
	echo $'\n\n'"Processing SV calling in chromosome $ii"
    Rscript ./scripts_to_publish/Clean_Manta.R $STRAIN $ii 1000000 0 > ./seqs/5_Clean_Manta/${ii}_log.txt
done