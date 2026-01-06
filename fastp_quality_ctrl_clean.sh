#!/bin/bash

# Fastp of Illumina reads with the ID of a strain
# Input arguments reads forward and reverse in fastq.gz
STRAIN=$1
SEQ_NGS_R1=$2
SEQ_NGS_R2=$3

mkdir -p seqs/1_qual_control

fastp -i "$SEQ_NGS_R1" -I "$SEQ_NGS_R2" -o ./seqs/1_qual_control/${STRAIN}_R1_fp_clean.fastq.gz -O ./seqs/1_qual_control/${STRAIN}_R2_fp_clean.fastq.gz -h ./seqs/1_qual_control/report_fastp_${STRAIN}.html -w 8