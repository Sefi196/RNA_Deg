#!/bin/bash

#Files
dir="path_to_bams_aligned_to_transcriptome" 
genome="path_to_bams_aligned_to_genome"
transcripts.fa="path_to_transcripts.fa_file" 
gencode.v31.annotation.gtf="path-to_gencode.v31.annotation.gtf_file"

#Command used for quatifaction with Nanocount. 

for file in "$dir"/*.bam; do
  NanoCount -i ${file} --extra_tx_info -o "${file##*/}.tsv" -b "${file##*/}.bam" 
done

#Command used for quatifaction with Salmon. 

for file in "$dir"/*.bam; do
  salmon quant -p 8 -t ${transcripts.fa} -l A -a ${file} -o ${file##*/} --noEffectiveLengthCorrection
done

#Command used for quatifaction with feature counts. 

/subread/bin/featureCounts -L -a ${gencode.v31.annotation.gtf} -o /feature_counts/qaunt.ft ${genome} --primary
  
