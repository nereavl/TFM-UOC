#!/bin/bash

# Index and align reads using bwa-mem2 with the D. subobscura hologenome
HOLOGENOME=$1
STRAIN=$2 # Strain name
SEQ_NGS_R1=$3
SEQ_NGS_R2=$4
LIBRARY=$5 # Library preparation used for sequencing
PLUNIT=$6 # Platformo Unit

mkdir -p seqs/2_bwa_mem_align

# Index the hologenome
echo "Indexing the hologenome with bwa-mem2..."
bwa-mem2 index $HOLOGENOME

# Align (map) reads to the hologenome and generate a SAM file
echo "Aligning reads to the hologenome..."
bwa-mem2 mem -t 8 -M -R "@RG\tID:D1EFBACXX_${STRAIN}\tPL:ILLUMINA\tLB:${LIBRARY}\tPU:${PLUNIT}\tSM:${STRAIN}" "$HOLOGENOME" "$SEQ_NGS_R1" "$SEQ_NGS_R2" > ./seqs/2_bwa_mem_align/${STRAIN}_aligned_hologenome_dsub.sam 2> ./seqs/2_bwa_mem_align/bwa_log.txt

