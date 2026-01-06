#!/bin/bash

# Processing SAM/BAM files from the mapping (alignment) of OF74 reads using samtools
STRAIN=$1
SAM_ALIGNMENT=$2

mkdir -p seqs/3_BAM_file

echo "Converting SAM to BAM..."
samtools view -@ 10 -b "$SAM_ALIGNMENT" -o "./seqs/3_BAM_file/${STRAIN}_mapped_bwamem2.bam"

echo "Grouping reads by name..."
samtools collate -@ 10 "./seqs/3_BAM_file/${STRAIN}_mapped_bwamem2.bam" -o "./seqs/3_BAM_file/${STRAIN}_collated.bam" 

echo "Fixing mate pair information..."
samtools fixmate -m "./seqs/3_BAM_file/${STRAIN}_collated.bam" "./seqs/3_BAM_file/${STRAIN}_fixmate.bam"

echo "Sorting BAM by coordinates..."
samtools sort -@ 10 -m 2G "./seqs/3_BAM_file/${STRAIN}_fixmate.bam" -o "./seqs/3_BAM_file/${STRAIN}_sorted.bam" 

echo "Marking and removing duplicates..."
samtools markdup -@ 10 "./seqs/3_BAM_file/${STRAIN}_sorted.bam" "./seqs/3_BAM_file/${STRAIN}_markdup_bwamem2.bam"
samtools markdup -@ 10 -r "./seqs/3_BAM_file/${STRAIN}_sorted.bam" "./seqs/3_BAM_file/${STRAIN}_remdup_bwamem2.bam"

echo "Indexing BAM file without duplicates..."
samtools index "./seqs/3_BAM_file/${STRAIN}_markdup_bwamem2.bam"
samtools index "./seqs/3_BAM_file/${STRAIN}_remdup_bwamem2.bam"

# Remove intermediate files to save storage space
rm -f "./seqs/3_BAM_file/${STRAIN}_mapped_bwamem2.bam" "./seqs/3_BAM_file/${STRAIN}_collated.bam" "./seqs/3_BAM_file/${STRAIN}_fixmate.bam" "./seqs/3_BAM_file/${STRAIN}_sorted.bam"

# Generate final BAM statistics
echo "Generating statistics for the final BAM file..."
samtools flagstat "./seqs/3_BAM_file/${STRAIN}_markdup_bwamem2.bam" > "./seqs/3_BAM_file/${STRAIN}_stats_markdup.txt"
samtools flagstat "./seqs/3_BAM_file/${STRAIN}_remdup_bwamem2.bam" > "./seqs/3_BAM_file/${STRAIN}_stats_remdup.txt"
