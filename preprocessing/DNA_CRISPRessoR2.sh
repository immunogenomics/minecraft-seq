#!/bin/bash
# Date: 02/06/2025
# Author: Michelle Curtis, Yuriy Baglaenko
# This script performs alignment of DNA FASTQs to an amplicon of interest using CRISPResso. 

# sed -i 's/\r//' Plate_amplicons.txt

OUT_DIR="/data/srlab1/mcurtis/GSK/Experiments/20240924_FBXO11/DNAOutput_R2/"
BATCH_BASE="/data/srlab1/mcurtis/GSK/Experiments/CRISPRessoBatch_R2.txt" 

while read -r PLATE SEQ;
do
    if [ -f "${OUT_DIR}${PLATE}/CRISPRessoBatch_on_batch/Nucleotide_percentage_summary.txt" ];    # Skip any plates already done
    then
        echo "Skipping $PLATE"
        continue
    fi

    FASTQ_DIR="/data/srlab1/GSK/Data/FBXO11/1726503765/fastq/Demultiplex/${PLATE}/"
    OUTPUT="${OUT_DIR}${PLATE}"
    mkdir -p $OUTPUT
    cd $OUTPUT

    PLATE_BATCH="${OUTPUT}/batch.batch"
    
    head -1 $BATCH_BASE > $PLATE_BATCH

    ### Subset batch files to include only the files present in directory
    ls $FASTQ_DIR*_R2.fastq.gz >> dna_filepaths.txt
    sed '2,$s/.fastq/.fastq.gz/g' $BATCH_BASE | awk -v FASTQ_DIR="$FASTQ_DIR" -v OFS="\t" '$2=FASTQ_DIR$2' | awk 'NR!=1 {print}' | grep -f dna_filepaths.txt >> $PLATE_BATCH
    rm dna_filepaths.txt

    CRISPRessoBatch -o ${OUTPUT} --batch_settings ${PLATE_BATCH}  --skip_failed -p 4 -a ${SEQ} --bo ${OUTPUT}  

done < Plate_amplicons.txt
