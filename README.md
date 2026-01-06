## **Requirements**
To run this pipeline, ensure the following tools are installed:

- **`fastp`**: Quality control and filtering of sequencing reads.  
  - Installation: [https://github.com/OpenGene/fastp](https://github.com/OpenGene/fastp)  

- **`bwa-mem2`**: Mapping reads to the reference genome.  
  - Installation: [https://github.com/bwa-mem2/bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)   

- **`Samtools`**: Processing and manipulation of SAM/BAM files.  
  - Installation: [https://www.htslib.org/doc/samtools.html](https://www.htslib.org/doc/samtools.html)

- **`Breakdancer`**: Detection of structural variants and breakpoints.  
  - Installation: [https://github.com/genome/breakdancer](https://github.com/genome/breakdancer)
 
- **`Lumpy`**: Detection of structural variants and breakpoints.  
  - Installation: 

- **`Manta`**: Detection of structural variants and breakpoints.  
  - Installation:
 
- **`GRIDSS`**: Detection of structural variants and breakpoints.  
  - Installation:

- **`R`**: Data analysis and filtering of structural variants
  - Installation:
  - Packages:
    - **`dplyr`**:
	- **`StructuralVariantAnnotation`**:
	

## **Workflow Summary**

Preparation
0. **Hologenome generation**: This script generates a hologenome by concatenating genomes from different species. It also processes the reference genome of *Drosophila subobscura* to ensure compatibility with downstream analyses.

   - **Script:** [`Hologenome_RefGenome_Dsubobscura.sh`](Hologenome_RefGenome_Dsubobscura.sh)  
   - **Command:**
	```bash
	./Create_Dsub_hologenome/Hologenome_RefGenome_Dsubobscura.sh
	```

The main pipeline follows these main steps:

1. **Quality Control**: Filtering and trimming of raw Illumina reads using `fastp`.  
   - **Script:** [`fastp_quality_ctrl_clean.sh`](fastp_quality_ctrl_clean.sh)  
   - **Command:**
     ```bash
     ./Scripts_workflow/fastp_quality_ctrl_clean.sh <STRAIN> <R1> <R2>
     ```
**Required arguments**
```
     -STRAIN  name of the strain of study (string)
     -R1  file containing R1 raw sequences (fastq.gz)
     -R2  file containing R2 raw sequences (fastq.gz)
```

2. **Read Mapping**: Aligning reads to the *Drosophila subobscura* hologenome using `bwa-mem2`.  
   - **Script:** [`Map_bwamem2_hologenome_Dsub.sh`](Map_bwamem2_hologenome_Dsub.sh)  
   - **Command:**
     ```bash
	 ./Scripts_workflow/Map_bwamem2_hologenome_Dsub.sh <ref_genome> <STRAIN> <R1> <R2> <LIBRARY> <PLATFORM_UNIT>
     ```
**Required arguments**
```
     -ref_genome  reference genome file (.fa)
     -STRAIN  name of the strain of study (string)
     -R1  file containing R1 filtered sequences (fastq.gz)
     -R2  file containing R2 filtered sequences (fastq.gz)
	 -LIBRARY  name of the library (string)
	 -PLATFORM_UNIT  name of the platform unit (string)
```

3. **BAM Processing**: Sorting, duplicate removal, and indexing using `samtools`.  
   - **Script:** [`Process_SAM_BAM_dup_removed.sh`](Process_SAM_BAM_dup_removed.sh)  
   - **Command:**
     ```bash
     ./Process_SAM_BAM_dup_removed.sh <STRAIN> <SAM_FILE>
     ```
**Required arguments**
	```
     -STRAIN  name of the strain of study (string)
     -SAM_FILE  raw SAM file (.sam)
	```

4.1 **Structural Variant Detection with BREAKDANCER**: This script prepares files for the identification of breakpoints and structural variants using `Breakdancer`. It also filters the results. 
   - **Script:** [`Process_Breakdancer.sh`](Process_Breakdancer.sh)  
   - **Command:**
     ```bash
     Process_Breakdancer.sh [-y] [-q] <STRAIN> <BAM_FILE> <CHR_LIST> <min_length> <min_quality>
     ```
**Requirements**
	```
     -BreakDancer
     -samtools
     -R (Clean_Breakdancer.R script)
	```
**Required arguments**
	```
     -STRAIN  name of the strain of study (string)
     -BAM_FILE  processed bam file (.bam)
     -CHR_LIST  list of chromosome names (string)
	```

4.2 **Structural Variant Detection with LUMPY**: Identification of breakpoints and structural variants using `LUMPY`.  
   - **Script:** [`Process_Lumpy.sh`](Process_Lumpy.sh)  
   - **Command:**
     ```bash
     Process_Lumpy.sh <STRAIN> <BAM_FILE> <ref_genome> <CHR_LIST>
     ```
**Required arguments**
	```
     -STRAIN  name of the strain of study (string)
     -BAM_FILE  processed bam file (.bam)
     -ref_genome  reference genome file (.fa)
     -CHR_LIST  list of chromosome names (string)
	```

4.3 **Structural Variant Detection with Manta**: Identification of breakpoints and structural variants using `Manta`.  
   - **Script:** [`Process_Manta.sh`](Process_Breakdancer.sh)  
   - **Command:**
     ```bash
     Process_Manta.sh <STRAIN> <BAM_FILE> <ref_genome> <CHR_LIST>
     ```
**Required arguments**
	```
     -STRAIN  name of the strain of study (string)
     -BAM_FILE  processed bam file (.bam)
     -ref_genome  reference genome file (.fa)
     -CHR_LIST  list of chromosome names (string)
	```
	
4.4 **Structural Variant Detection with GRIDSS**: Identification of breakpoints and structural variants using `GRIDSS`.  
   - **Script:** [`Process_GRIDSS.sh`](Process_GRIDSS.sh)  
   - **Command:**
     ```bash
     Process_GRIDSS.sh <STRAIN> <BAM_FILE> <ref_genome> <CHR_LIST>
     ```
**Required arguments**
	```
     -STRAIN  name of the strain of study (string)
     -BAM_FILE  processed bam file (.bam)
     -ref_genome  reference genome file (.fa)
     -CHR_LIST  list of chromosome names (string)
	```		
	
Each step generates key output files that contribute to detecting chromosomal inversions in *Drosophila subobscura*.





