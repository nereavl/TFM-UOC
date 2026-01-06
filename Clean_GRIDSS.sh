#!/bin/bash

# Clean the outputs by breakdancer
# The generated txt file is stored in a separate folder to improve result visibility.

STRAIN=$1
DSUB_CHR=$2  # RefSeq chromosome code (e.g., NC_048534.1)
MIN_LENTGH=$3
MIN_QUALITY=$4

mkdir -p ./seqs/5_Clean_GRIDSS

for ii in $DSUB_CHR; do
	echo $'\n\n'"Processing SV calling in chromosome $ii"
    Rscript ./scripts/Scripts_workflow/Clean_GRIDSS.R $STRAIN $ii $MIN_LENTGH $MIN_QUALITY > ./seqs/5_Clean_GRIDSS/${ii}_log.txt
done
