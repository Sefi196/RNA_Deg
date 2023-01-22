#!/bin/bash
#SBATCH --partition="physical"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128000
#time in days-hours:mins:sec
#SBATCH --time=1-00:00:00

# Takes roughly 1hr for 5mil reads

dir="all_fastq"

for f in "dir"
do
echo "Processing $f"
./scripts/dRNA_mapping_quant hg38.fasta gencode.v31.transcripts.fa gencode.v31.annotation.gtf sequin_genome.fasta sequin_transcriptome.fasta sequin_annotation.gtf path/to/fastq/$f/pass.fastq $f

done
