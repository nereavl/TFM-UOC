#!/bin/bash

# Input as an agument: tab-delimited file containing the windows of the regions of interest and cromosome (RefSeq: NC_048534.1, for ex), sequences
# FASTQ.gz cleaned, the BAM alignment and the reference genome

TAB_FILE_IN=$1
SEQ_NGS_R1=$2
SEQ_NGS_R2=$3
BAM_ALIGNMENT=$4
REF_GENOME=$5 # Reference genome in FASTA format

# IFS=$'\t' to read the tab-delimited file line by line and create directories for reads and assembly (to obtain the contigs) for each window
while IFS=$'\t' read -r ID_window output_folder ID_chr proximal_start proximal_end distal_start distal_end; do

# The output folders for reads and assembly are created to keep the files organized and make it easier to manage the results.
 echo $"Create Folders"
 reads_folder="./seqs/6_denovo_assembly/${output_folder}/reads_split"
 assembly_folder="./seqs/6_denovo_assembly/${output_folder}/assembly_contigs"
 mkdir -p "./seqs/6_denovo_assembly/${output_folder}" "$reads_folder" "$assembly_folder"

# Extract reads for the chromosome U breakpoint regions (RefSeq: NC_048534.1). Variable names indicate whether they correspond to the proximal or distal breakpoint.
# -F 12 ensures that only properly mapped reads at the breakpoints are extracted.
 echo $"Extract reads from breakpoint regions"
 samtools view -F 12 -h -b "$BAM_ALIGNMENT" "Drosophila_subobscura_${ID_chr}:${proximal_start}-${proximal_end}" > "${reads_folder}/${ID_window}_mapped_proximal.bam"
 samtools view -F 12 -h -b "$BAM_ALIGNMENT" "Drosophila_subobscura_${ID_chr}:${distal_start}-${distal_end}" > "${reads_folder}/${ID_window}_mapped_distal.bam"

# Extract discordant reads where one end is mapped, but the other is not
 echo $"Extract discordant reads where one end is mapped, but the other is not"
 samtools view -F 4 -f 8 -b "$BAM_ALIGNMENT" "Drosophila_subobscura_${ID_chr}:${proximal_start}-${proximal_end}" > "${reads_folder}/${ID_window}_discordant_proximal_1.bam"
 samtools view -F 8 -f 4 -b "$BAM_ALIGNMENT" "Drosophila_subobscura_${ID_chr}:${proximal_start}-${proximal_end}" > "${reads_folder}/${ID_window}_discordant_proximal_2.bam"
 samtools view -F 4 -f 8 -b "$BAM_ALIGNMENT" "Drosophila_subobscura_${ID_chr}:${distal_start}-${distal_end}" > "${reads_folder}/${ID_window}_discordant_distal_1.bam"
 samtools view -F 8 -f 4 -b "$BAM_ALIGNMENT" "Drosophila_subobscura_${ID_chr}:${distal_start}-${distal_end}" > "${reads_folder}/${ID_window}_discordant_distal_2.bam"

# Merge reads from all SAM files, including properly mapped and discordant reads
 echo $"Merge reads from all SAM files, including properly mapped and discordant reads"
 samtools merge -f "${reads_folder}/${ID_window}_final.bam" \
 "${reads_folder}/${ID_window}_mapped_proximal.bam" "${reads_folder}/${ID_window}_mapped_distal.bam" \
 "${reads_folder}/${ID_window}_discordant_proximal_1.bam" "${reads_folder}/${ID_window}_discordant_proximal_2.bam" \
 "${reads_folder}/${ID_window}_discordant_distal_1.bam" "${reads_folder}/${ID_window}_discordant_distal_2.bam"

# Sort the BAM file by read name
 echo $"Sort the BAM file by read name"
 samtools sort -n -o "${reads_folder}/${ID_window}_final_sorted.bam" "${reads_folder}/${ID_window}_final.bam"

# Split the reads into FASTQ files to extract singletons and get the unmapped reads from the original fastq files (cleaned)
 echo $"Split the reads into FASTQ files to extract singletons and get the unmapped reads from the original fastq files (cleaned)"
 samtools fastq -1 "${reads_folder}/${ID_window}_R1.fastq" -2 "${reads_folder}/${ID_window}_R2.fastq" \
 -s "${reads_folder}/${ID_window}_singles.fastq" "${reads_folder}/${ID_window}_final_sorted.bam"

# Extract the headers from singletons
 echo $"Extract the headers from singletons"
 grep '^@HWI' "${reads_folder}/${ID_window}_singles.fastq" > "${reads_folder}/names_singles.txt"

# Filter the readers that match with the extracted headers to get those reads that are not mapped with bwa-mem
 echo $"Filter the readers that match with the extracted headers to get those reads that are not mapped with bwa-mem"
 zcat "$SEQ_NGS_R1" | grep -A3 -wFf <(sed 's/\/[12]$//' "${reads_folder}/names_singles.txt") > "${reads_folder}/unmapped_R1.fastq"
 zcat "$SEQ_NGS_R2" | grep -A3 -wFf <(sed 's/\/[12]$//' "${reads_folder}/names_singles.txt") > "${reads_folder}/unmapped_R2.fastq"

 sed '/^--/d' "${reads_folder}/unmapped_R1.fastq" > "${reads_folder}/clean_unmapped_R1.fastq"
 sed '/^--/d' "${reads_folder}/unmapped_R2.fastq" > "${reads_folder}/clean_unmapped_R2.fastq"

# Combine retrieved reads into final input files R1 and R2 for SPAdes
 echo $"Combine retrieved reads into final input files R1 and R2 for SPAdes"
 cat "${reads_folder}/${ID_window}_R1.fastq" "${reads_folder}/clean_unmapped_R1.fastq" > "${reads_folder}/${ID_window}_final_R1.fastq"
 cat "${reads_folder}/${ID_window}_R2.fastq" "${reads_folder}/clean_unmapped_R2.fastq" > "${reads_folder}/${ID_window}_final_R2.fastq"

# Run SPAdes assembly and save the results in a separated folder
 echo $"Run SPAdes assembly and save the results in a separated folder"
 spades.py -1 "${reads_folder}/${ID_window}_final_R1.fastq" -2 "${reads_folder}/${ID_window}_final_R2.fastq" -o "$assembly_folder" --careful

# Contigs can be filtered with BLAST using the reference genome passed as an argument "$REF_GENOME"
 echo $"BLAST contigs against reference genome"
 blastn -query "${assembly_folder}/contigs.fasta" -subject "$REF_GENOME" -out "./seqs/6_denovo_assembly/${output_folder}/${ID_window}_blast_results.txt" -outfmt 6

# From previous BLASTn web tests, contigs sometimes align to other chromosomes with low similarity.
# Filter BLAST results to  include matches to chromosome selected in the tab file
 echo $"Filter BLAST results"
 grep "Drosophila_subobscura_${ID_chr}" "./seqs/6_denovo_assembly/${output_folder}/${ID_window}_blast_results.txt" > "./seqs/6_denovo_assembly/${output_folder}/${ID_window}_blast_filtered.txt"

# Extend the flanking regions by 5Kb before and after the studied regions to ensure the full region of interest is captured
 awk '($9 >= ('$proximal_start' - 5000) && $9 <= ('$proximal_end' + 5000)) || ($9 >= ('$distal_start' - 5000) && $9 <= ('$distal_end' + 5000))' \
 "./seqs/6_denovo_assembly/${output_folder}/${ID_window}_blast_filtered.txt" > "./seqs/6_denovo_assembly/${output_folder}/${ID_window}_blast_flanking.txt"

done < "$TAB_FILE_IN"
