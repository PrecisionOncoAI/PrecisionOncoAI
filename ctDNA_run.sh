#!/bin/bash

# Define input arguments
FASTQ1=$1
FASTQ2=$2
REFERENCE=$3
OUTPUT_DIR=$4

mkdir -p ${OUTPUT_DIR}

# Step 1: Quality Control
fastqc ${FASTQ1} ${FASTQ2} -o ${OUTPUT_DIR}/qc/

# Step 2: Trimming
trimmomatic PE -phred33 ${FASTQ1} ${FASTQ2} \
    ${OUTPUT_DIR}/trimmed_1.fastq.gz ${OUTPUT_DIR}/unpaired_1.fastq.gz \
    ${OUTPUT_DIR}/trimmed_2.fastq.gz ${OUTPUT_DIR}/unpaired_2.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# Step 3: Alignment
bwa mem -t 8 ${REFERENCE} ${OUTPUT_DIR}/trimmed_1.fastq.gz ${OUTPUT_DIR}/trimmed_2.fastq.gz > ${OUTPUT_DIR}/aligned.sam
samtools view -bS ${OUTPUT_DIR}/aligned.sam > ${OUTPUT_DIR}/aligned.bam
samtools sort ${OUTPUT_DIR}/aligned.bam -o ${OUTPUT_DIR}/sorted.bam
samtools index ${OUTPUT_DIR}/sorted.bam

# Step 4: Variant Calling
gatk HaplotypeCaller -R ${REFERENCE} -I ${OUTPUT_DIR}/sorted.bam -O ${OUTPUT_DIR}/variants.vcf

# Step 5: Annotation (Using a Python script)
python3 annotate_variants.py ${OUTPUT_DIR}/variants.vcf ${OUTPUT_DIR}/annotated_variants.txt
