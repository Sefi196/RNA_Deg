library(dplyr)
library(GenomicFeatures)
library(pheno2geno)
library(devtools)

#to filter gtf file for single isoform genes and to get 
#1. list of singe isofrom genes and trabscript names for filtering BamSLAM output. 
#2. bed file with gene coordiantes 

# read in gtf and pull length data
txs <- makeTxDbFromGFF("gencode.v31.annotation.gtf", format="gtf")
txLengths <- transcriptLengths(txs, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
lengths <- data.frame(txLengths$tx_name, txLengths$tx_len)

#filter GTF for single isofrm genes
single_isoform <- txLengths %>%
  group_by(gene_id) %>%
  dplyr::filter(length(unique(tx_name)) == 1) %>%
  dplyr::select(2,3)

# filter bamSLAM outputs 
single_isoform_gene_plotme <- merge(single_isoform, bamSLAM_output.csv, by.x="tx_name", by.y="transcript", all.x=FALSE)

#now can plot coverages and perform architecture analysis

########
#get bed file
#######

gff_file <- "pathto_gencode.v31.annotation.gtf"

gff_annot <- rtracklayer::import(gff_file) %>%
  as.data.frame() %>%
  tibble::as_tibble() 


#load gene file
gene_file <- singel_isofrom$gene_id


gff_filtered <- gff_annot %>%
  dplyr::filter(gene_id %in% genes_to_filter)

#need the (1,2,3) and the gene only

bed <- dplyr::filter(gff_filtered, , type == "gene") %>%
  dplyr::select(1,2,3)

write.table(bed, "single_isoform_position.bed", row.names = FALSE, quote = FALSE,sep="\t")
write.csv(single_isoform$tx_name, "single_isoform_genes_ENST.csv", row.names = FALSE, quote = FALSE)
write.table(singel_isofrom$gene_id, "single_isform_gene.list.csv", row.names = FALSE, quote = FALSE,sep="\t")
