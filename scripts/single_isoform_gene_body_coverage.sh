#!/bin/bash

#require samtools and bedtools
# samtools/1.11
# bedtools/2.30.0

dir="bams/" # sample target dir
bed="single_isoform_position.bed"
output="path_to_output_dir"


echo running loop
COUNTER=0

#run in loop on all bam files 
for file in "$dir"/*.bam; do
   #samtools view -h -b -L ${bed} ${file} > ${output}"filt_${file##*/}" #find overlapps
   bedtools intersect -abam ${file} -b ${bed} -f 1.0 -wa > ${output}"filt_${file##*/}" #find reads within
   samtools index ${output}"filt_${file##*/}"
   COUNTER=$((COUNTER+1))
   echo "Number of files run: $COUNTER"
done


bed_RSEQC="path_to_hg38_Gencode_V24.bed"
bams=""

geneBody_coverage.py -r ${bed_RSEQC} -i ${output} -o /data/gpfs/projects/punim1441/Project_dRNA_deg/single_isoform_bams/RSEQC/single_isoform_genes
