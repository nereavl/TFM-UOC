#!/bin/bash


############################################
# Script: run_breakdancer.sh
# Description: Detects chromosomal breakpoints
#              using LUMPY 
# Author: Mercé Merayo & Kenia Delgado
# Adaptation: Nerea Vallejo
############################################

# ----------------------------
# Positional arguments
# ----------------------------
STRAIN=$1
BAM_ALIGNMENT=$2
HOLOGENOME=$3
DSUB_CHR=$4

# ----------------------------
# Prepare folders
# ----------------------------
echo $"Preparing folders"
# mkdir -p ./seqs/3_SAM_for_Lumpy
# mkdir -p ./seqs/4_Lumpy
# mkdir -p ./seqs/5_Clean_Lumpy

# copy and index chosen BAM alignment
echo $"Copying and indexing BAM alignment"
# cp $BAM_ALIGNMENT ./seqs/3_SAM_for_Lumpy/${STRAIN}_mapped.bam
# samtools index ./seqs/3_SAM_for_Lumpy/${STRAIN}_mapped.bam

# copy and index chosen BAM hologenome file
echo $"Copying and indexing hologenome"
# cp $HOLOGENOME ./seqs/3_SAM_for_Lumpy/hologenome.fa
# samtools faidx ./seqs/3_SAM_for_Lumpy/hologenome.fa --write-index

# discordants
echo $"Extracting discordants"
# samtools view -@ 10 -b -F 1294 ./seqs/3_SAM_for_Lumpy/${STRAIN}_mapped.bam > ./seqs/3_SAM_for_Lumpy/${STRAIN}_discordants_unsorted.bam
# samtools sort -@ 10 -m 2G -o ./seqs/3_SAM_for_Lumpy/${STRAIN}_discordants.bam ./seqs/3_SAM_for_Lumpy/${STRAIN}_discordants_unsorted.bam
# samtools index ./seqs/3_SAM_for_Lumpy/${STRAIN}_discordants.bam
# rm ./seqs/3_SAM_for_Lumpy/${STRAIN}_discordants_unsorted.bam

# splitters
echo $"Extracting splitters"
# samtools view -h ./seqs/3_SAM_for_Lumpy/${STRAIN}_mapped.bam | scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ./seqs/3_SAM_for_Lumpy/${STRAIN}_splitters_unsorted.bam
# samtools sort -@ 10 -m 2G -o ./seqs/3_SAM_for_Lumpy/${STRAIN}_splitters.bam ./seqs/3_SAM_for_Lumpy/${STRAIN}_splitters_unsorted.bam
# samtools index ./seqs/3_SAM_for_Lumpy/${STRAIN}_splitters.bam
# rm ./seqs/3_SAM_for_Lumpy/${STRAIN}_splitters_unsorted.bam

# samtools view -r readgroup1 ./seqs/3_SAM_for_Lumpy/${STRAIN}_mapped.bam \
#     | tail -n+1000000 \
#     | scripts/pairend_distro.py \
#     -r 101 \
#     -X 4 \
#    -N 10000 \
#     -o ./seqs/3_SAM_for_Lumpy/OF74.lib1.histo


# Making BED files to exclude SV calling regions

for ii in $DSUB_CHR; do

  echo $"Preparing BED files for region exclusion for chromosome $ii"
	awk '{print $1, 0, $2}' OFS="\t" ./seqs/3_SAM_for_Lumpy/hologenome.fa.fai > ./seqs/3_SAM_for_Lumpy/full.bed
	grep -v "Drosophila_subobscura_$ii" ./seqs/3_SAM_for_Lumpy/full.bed > ./seqs/3_SAM_for_Lumpy/$ii.bed
	
	
	echo $'\n\n'"Detecting breakpoints in chromosome $ii"
	echo $"Running Lumpy"
	lumpyexpress -B ./seqs/3_SAM_for_Lumpy/${STRAIN}_mapped.bam -S ./seqs/3_SAM_for_Lumpy/${STRAIN}_splitters.bam -D ./seqs/3_SAM_for_Lumpy/${STRAIN}_discordants.bam -x ./seqs/3_SAM_for_Lumpy/$ii.bed -o ./seqs/4_Lumpy/${STRAIN}_lumpy_$ii.vcf

#     	lumpy \
#     		-mw 4 \
#     		-tt 0 \
#     		-x ./seqs/3_SAM_for_Lumpy/$ii.bed \
#     		-pe id:${ii},bam_file:./seqs/3_SAM_for_Lumpy/${STRAIN}_discordants.bam,histo_file:./seqs/3_SAM_for_Lumpy/${STRAIN}.lib1.histo,mean:306,stdev:67,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:0 \
#     		-sr id:sample,bam_file:./seqs/3_SAM_for_Lumpy/${STRAIN}_splitters.bam,back_distance:10,weight:1,min_mapping_threshold:0 \
#     		> ./seqs/4_Lumpy/${STRAIN}_lumpy_$ii.vcf
	
	echo $'\n\n'"Processing SV calling in chromosome $ii"
	Rscript ./scripts_to_publish/Clean_LUMPY.R $STRAIN $ii 1000000 0 > ./seqs/5_Clean_Lumpy/${ii}_log.txt
	
done
