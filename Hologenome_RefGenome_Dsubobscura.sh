#!/bin/bash

# List of symbiotic and pathogenic species of D. subobscura, including this species itself
# with the download location URL-FTP of the NCBI reference genomes

list_species=(
# NCBI reference genomes
"Saccharomyces_cerevisiae|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"
"Wolbachia_pipientis|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/971/685/GCF_007971685.1_ASM797168v1/GCF_007971685.1_ASM797168v1_genomic.fna.gz"
"Pseudomonas_entomophila|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/105/GCF_000026105.1_ASM2610v1/GCF_000026105.1_ASM2610v1_genomic.fna.gz"
"Commensalibacter_intestini|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/231/445/GCF_000231445.1_ASM23144v1/GCF_000231445.1_ASM23144v1_genomic.fna.gz"
"Acetobacter_pomorum|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/738/225/GCF_002738225.1_ASM273822v1/GCF_002738225.1_ASM273822v1_genomic.fna.gz"
"Gluconobacter_morbifer|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/234/355/GCF_000234355.1_ASM23435v1/GCF_000234355.1_ASM23435v1_genomic.fna.gz"
"Providencia_burhodogranariea|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/314/855/GCF_000314855.2_ASM31485v2/GCF_000314855.2_ASM31485v2_genomic.fna.gz"
"Providencia_alcalifaciens|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/393/505/GCF_002393505.1_ASM239350v1/GCF_002393505.1_ASM239350v1_genomic.fna.gz"
"Providencia_rettgeri|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/048/105/GCF_019048105.1_ASM1904810v1/GCF_019048105.1_ASM1904810v1_genomic.fna.gz"
"Enterococcus_faecalis|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/393/015/GCF_000393015.1_Ente_faec_T5_V1/GCF_000393015.1_Ente_faec_T5_V1_genomic.fna.gz"
"Lactiplantibacillus_plantarum|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/913/655/GCF_009913655.1_ASM991365v1/GCF_009913655.1_ASM991365v1_genomic.fna.gz"
"Levilactobacillus_brevis|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/475/625/GCF_900475625.1_45709_F01/GCF_900475625.1_45709_F01_genomic.fna.gz"
"Hanseniaspora_uvarum|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/037/102/615/GCA_037102615.1_ASM3710261v1/GCA_037102615.1_ASM3710261v1_genomic.fna.gz"
"Pichia_kudriavzevii|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/054/445/GCF_003054445.1_ASM305444v1/GCF_003054445.1_ASM305444v1_genomic.fna.gz"
"Maudiozyma_humilis|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/037/102/785/GCA_037102785.1_ASM3710278v1/GCA_037102785.1_ASM3710278v1_genomic.fna.gz"
"Candida_sake|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/243/815/GCA_003243815.1_ASM324381v1/GCA_003243815.1_ASM324381v1_genomic.fna.gz"
"Rhodotorula_babjevae|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/920/104/245/GCA_920104245.2_Hybrid_genome_assembly_and_annotation/GCA_920104245.2_Hybrid_genome_assembly_and_annotation_genomic.fna.gz"
"Serratia_liquefaciens|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/422/085/GCF_000422085.1_ASM42208v1/GCF_000422085.1_ASM42208v1_genomic.fna.gz"
"Trigonopsis_vinaria|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/582/095/GCA_030582095.1_ASM3058209v1/GCA_030582095.1_ASM3058209v1_genomic.fna.gz"
"Acetobacter_aceti|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/379/545/GCF_000379545.1_ASM37954v1/GCF_000379545.1_ASM37954v1_genomic.fna.gz"
"Leuconostoc_mesenteroides|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/014/445/GCF_000014445.1_ASM1444v1/GCF_000014445.1_ASM1444v1_genomic.fna.gz"
"Enterobacter_cloacae|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/905/331/265/GCF_905331265.2_AI2999v1_cpp/GCF_905331265.2_AI2999v1_cpp_genomic.fna.gz"
"Klebsiella_oxytoca|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/636/985/GCF_900636985.1_45889_C01/GCF_900636985.1_45889_C01_genomic.fna.gz"
"Pectobacterium_carotovorum|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/488/025/GCF_013488025.1_ASM1348802v1/GCF_013488025.1_ASM1348802v1_genomic.fna.gz"
"Zygosaccharomyces_bailii|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/949/129/075/GCA_949129075.1_ZBA_7846_DN/GCA_949129075.1_ZBA_7846_DN_genomic.fna.gz"

# Reference genome of Drosophila subobscura
"Drosophila_subobscura|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/121/235/GCF_008121235.1_UCBerk_Dsub_1.0/GCF_008121235.1_UCBerk_Dsub_1.0_genomic.fna.gz"
)

mkdir ref_genome

for each_specie in "${list_species[@]}"; do
 specie=$(echo "$each_specie" | cut -d'|' -f1)
 url=$(echo "$each_specie" | cut -d'|' -f2)
 curl -O "$url"
 file_specie_dwnld=$(basename "$url")
# Decompress and extract the .fasta file from the downloaded .gz file.
# In the new header (starting with >), the species name from the list is added
# along with the chromosome identifier (which is in the original header until the space) using `sed`
 gunzip -c "$file_specie_dwnld" | sed "s/^>/>${specie}_/; s/ .*$//" > "${specie}.fasta"
done

# Merge all reference genomes into a multi-FASTA file, resulting in a hologenome.
# The reference genome of *D. subobscura* is saved with a .fa extension to avoid deletion at the end of the script.
cat *.fasta > ref_genome/hologenome_with_dsub.fa
cat Drosophila_subobscura.fasta > ref_genome/Drosophila_subobscura_ref_genome.fa

# Delete the downloaded species files to free up memory. 
rm *.gz
rm *.fasta
