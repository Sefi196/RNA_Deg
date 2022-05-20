main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  countdata <- args[1]
  output <- args[2]
  suppressPackageStartupMessages({
    library(dplyr, warn.conflicts = FALSE)
    library(ggplot2)
    library(reshape)
    library(tidyr)
    library(DropletUtils)
    library(Hmisc)
  })
  
  #Import data 
  #Genes from NanoCount
  setwd("~/Documents/Ph.D/RNA_degradation_study/Deseq_15_samples/20210729_ft_counts/29_July_feature_counts/")
  CountMatrix <- read.csv("ft.counts.csv", row.names = 1, header = T)
  
  #GIsofroms from Feature counts 
  setwd("~/Documents/Ph.D/RNA_degradation_study/Deseq_15_samples/NanoCount-0.4.0/")
  CountMatrix <- read.csv("NanoCount.Count.Matrix.csv", row.names = 1, header = T)

  #CountMatrix <- read.csv(countdata, row.names = 1, header = T)
  message("Importing Count Matrix")
  CountMatrix <- as.matrix(as.data.frame(CountMatrix))
  colnames(CountMatrix) <- c("TS12_RIN_9.9", "TS10_RIN_9.8", "TS11_RIN_9.7", "TS10_RIN_9.6","TS11_RIN_9.6", "TS10_RIN_9.3", "TS11_RIN_9.3", "TS11_RIN_8.9", "TS11_RIN_8.8", 
                             "TS10_RIN_8.4", "TS11_RIN_8.7", "TS12_RIN_8.2", "TS10_RIN_7.7", "TS12_RIN_7.3", "TS12_RIN_7.2"
                             #, "TS12_RIN_9.9_cDNA", "TS12_RIN_9.9_2_cDNA", "TS11_RIN_9.8_cDNA", "TS12_RIN_8.2_cDNA", "TS12_RIN_8.2_2_cDNA", "TS10_RIN_7.7_cDNA"
                             )

  #remove rows that sum to zero
  CountMatrix <- CountMatrix[ rowSums(CountMatrix[,-1]) > 0, ]

  ################################################################################
  ## Subsampling - saturation
  ################################################################################
  
  message("Defining Functions")
  
  #Define Functions
  #Define a function that will add replicates together based on there Time to degradation treatment. 
  add_sum_time <-function(CountMatrix) {
   tmp <- data.frame(CountMatrix) %>%
      tibble::rownames_to_column("feature") %>%
      tidyr::gather(key = "sample", value = "count", -feature)  %>%
      dplyr::mutate("type" = 
                  case_when(sample == "TS12_RIN_9.9"  ~ 'CONTROL', 
                            sample == "TS10_RIN_9.8" ~ 'CONTROL',
                            sample == "TS11_RIN_9.7" ~ 'CONTROL',
                            sample == "TS10_RIN_9.6" ~ 'VERY HIGH',
                            sample == "TS11_RIN_9.6" ~ 'VERY HIGH',
                            sample == "TS10_RIN_9.3" ~ 'HIGH',
                            sample == "TS11_RIN_9.3" ~ 'HIGH',
                            sample == "TS11_RIN_8.9" ~ 'MID',
                            sample == "TS11_RIN_8.8" ~ 'MID',
                            sample == "TS10_RIN_8.4" ~ 'MID',
                            sample == "TS11_RIN_8.7" ~ 'LOW',
                            sample == "TS12_RIN_8.2" ~ 'LOW',
                            sample == "TS10_RIN_7.7" ~ 'LOW',
                            sample == "TS12_RIN_7.3" ~ 'VERY LOW',
                            sample == "TS12_RIN_7.2" ~ 'VERY LOW'
                            #sample == "TS12_RIN_9.9_cDNA"  ~ 'CONTROL cDNA',
                            #sample == "TS12_RIN_9.9_2_cDNA"  ~ 'CONTROL cDNA',
                            #sample == "TS11_RIN_9.8_cDNA"  ~ 'CONTROL cDNA',
                            #sample == "TS12_RIN_8.2_cDNA"  ~ 'LOW cDNA',
                            #sample == "TS12_RIN_8.2_2_cDNA"  ~ 'LOW cDNA',
                            #sample == "TS10_RIN_7.7_cDNA"  ~ 'LOW cDNA'
                            ))
    tmp %>%
      dplyr::group_by(type, feature) %>%
      dplyr::summarise(sum(count))
  } 


  #Preparing the data 
  #Use add sum function to merge replicates into one sample 
  temp <- add_sum_time(CountMatrix)

  #Define Function that will revert the merged counts into a count matrix for sub sampling and graphing 

  reverttoCountMatrix <-function(temp) {
    CountMatrix <- temp %>%
      tidyr::pivot_wider(names_from = "type", values_from = "sum(count)")
    CountMatrix <- data.frame(CountMatrix, row.names = 1)
    
  }
  
  message("Preparing data")

  #Preparing the data 
  #Use revert to count Matrix function 
  CountMatrix <- as.matrix(reverttoCountMatrix(temp))
  
  message("Define vector for subfractions")
  
  #define vector between 0 and 1 to be used for proportion of sub sampling of proportions 
  subfracs <- c(seq(from = 0.000, to = 0.025, length.out = 10), 
              c(0.035, 0.045, 0.06, 0.08), 
              seq(from = 0.1, to = 1, by = 0.05))


  #Function to sub sample count matrix 
  #pull out number of genes and number of reads and return a data frame. 
  subsample <- function(val) {
    genemat <- DropletUtils::downsampleMatrix(CountMatrix, prop = val, bycol = TRUE) 
    tmp123 <- data.frame(sample = colnames(CountMatrix),
             nGene = colSums(genemat >= 1),
             nReadsGene = colSums(genemat),
             stringsAsFactors = FALSE) %>%
   dplyr::mutate(sampleFrac = val)
  return(as.data.frame(tmp123))
  }

  #use lappy function to cycle through the subfracs and return data frame for every val (sub fraction proportion) 
  plot.me_subfracs <- as.data.frame(do.call(dplyr::bind_rows, lapply(subfracs, subsample)))
  
  message("Add time to lysis")

  #Define a function that will add Time as a col to my data frame 
  add_time <- function(df) {
    df %>% mutate(Time =
                      case_when(sample == "CONTROL"  ~ '0 hours', 
                                sample == "VERY.HIGH" ~ '0.5 hours',
                                sample == "HIGH" ~ '1 hours',
                                sample == "MID" ~ '3_4 hours',
                                sample == "LOW" ~ '6 hours',
                                sample == "VERY.LOW" ~ '8 hours',
                                sample == "CONTROL.cDNA"  ~ '0 hours',
                                sample == "LOW.cDNA" ~ '6 hours'
                                )
  )
  }
  
  #Define a function that will add method i.e dRNA vs cDNA
  add_method <- function(df) {
    df %>% mutate(Method =
                    case_when(grepl("cDNA", sample) ~ "cDNA",
                              !grepl("cDNA", sample) ~ "dRNA"
                              )) 
  }

  #Call the add_time function so we can plot with time as a group 
  plot.me_subfracs <- add_time(plot.me_subfracs) 
  plot.me_subfracs <- add_method(plot.me_subfracs) 

  #Set the levels 
  plot.me_subfracs$sample <- factor(plot.me_subfracs$sample, levels = c('CONTROL','VERY.HIGH','HIGH','MID','LOW','VERY.LOW'
                                                                        #, "CONTROL.cDNA", "LOW.cDNA"
                                                                        ))
  plot.me_subfracs$Time <- factor(plot.me_subfracs$Time, levels = c('0 hours','0.5 hours','1 hours','3_4 hours','6 hours','8 hours'))
  plot.me_subfracs$Method <- factor(plot.me_subfracs$Method, levels = c('dRNA'
                                                                        #, 'cDNA'
                                                                        ))


  #Use temp data to plot data saturation curves per sample 
  message("Plotting")
  
  plot.me_subfracs.dRNA <- filter(plot.me_subfracs, Method  == 'dRNA')
  
  pdf(paste0(output, "RIN_.pdf"), width=8, height=5)
  q1 <- ggplot(plot.me_subfracs,
             aes(x = nReadsGene, y = nGene, colour=sample, group = sample)) +
    geom_line(colour = "Black") + 
    geom_point(size=2.5) +
    theme_bw() + xlab("Sampled number of assigned reads") + 
    ggtitle("Genes") + 
    ylab("Number of detected features") +
    xlim(0,3000000) +
    #ylim(0,20000) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
    plot.title = element_text(color="black", size=20),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
    theme(legend.title = element_blank())
  print(q1)
  dev.off()
  
  #To plot control vs very.low 
  plot.me_subfracs.dRNA <- filter(plot.me_subfracs, sample  == 'CONTROL' | sample  == 'VERY.LOW')

  pdf(paste0(output, "time_.pdf"), width=8, height=5)
  pdf("Isoforms_Saturation_dRNA.pdf", width=8, height=5) # run this for now but need to get rid for clean function
  q2 <- ggplot(plot.me_subfracs.dRNA,
               aes(x = nReadsGene, y = nGene, colour=Time, group = Time)) +
    geom_line(colour = "Black") +
    geom_point(size=2.5) + 
    scale_colour_manual("legend", values = c("0 hours"="#F8766D", "8 hours"="#F564E3")) +
    theme_bw() + xlab("Sampled number of assigned reads") + 
    ggtitle("Genes") + 
    ylab("Number of detected features") +
    xlim(0,2000000) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
      plot.title = element_text(color="black", size=20),
      axis.title.x = element_text(color="black", size=14),
      axis.title.y = element_text(color="black", size=14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()) +
    theme(legend.title = element_blank())
  print(q2)
  dev.off()
  
  #plot cDNA vs dRNA
  #filter cDNA vs dRNA samples
  dRNAvcDNA <- filter(plot.me_subfracs, Method  == 'cDNA' | sample  == 'CONTROL' | sample  == 'LOW')
  cDNA <- filter(plot.me_subfracs, Method  == 'cDNA')
  
  
  pdf(paste0(output, "dRNAvcDNA_RIN_.pdf"), width=8, height=4)
  pdf("Saturation_dRNAvcDNA_RIN_.pdf", width=8, height=4) # run this for now but need to get rid for clean function
  q3 <- ggplot(dRNAvcDNA,
               aes(x = nReadsGene, y = nGene, colour=Method, group = Method)) +
    geom_line(colour = "Black") + 
    facet_grid(~ Time) +
    geom_point(size=2.5) +
    theme_bw() + xlab("Sampled number of assigned reads") + 
    ggtitle("Genes") + 
    ylab("Number of detected features") +
    xlim(0,2000000) +
    ylim(0,17000) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
      plot.title = element_text(color="black", size=20),
      axis.title.x = element_text(color="black", size=14),
      axis.title.y = element_text(color="black", size=14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()) +
    theme(legend.title = element_blank())
  print(q3)
  dev.off()
  
  
  message("Completed rareafaction curve")
  
  ################################################################################
  ## Detection by transcript length
  ################################################################################
  setwd("~/Documents/Ph.D/RNA_degradation_study/Deseq_15_samples/20210729_ft_counts/29_July_feature_counts/")
  CountMatrix <- read.csv("ft.counts.csv", row.names = 1, header = T)
  
  CountMatrix <- read.csv("PCS110_dRNA_ft.counts.csv", row.names = 1, header = T)
  
  CountMatrix <- as.matrix(as.data.frame(CountMatrix))
  colnames(CountMatrix) <- c("TS12_RIN_9.9", "TS10_RIN_9.8", "TS11_RIN_9.7", "TS10_RIN_9.6","TS11_RIN_9.6", "TS10_RIN_9.3", "TS11_RIN_9.3", "TS11_RIN_8.9", "TS11_RIN_8.8", 
                             "TS10_RIN_8.4", "TS11_RIN_8.7", "TS12_RIN_8.2", "TS10_RIN_7.7", "TS12_RIN_7.3", "TS12_RIN_7.2"
                             , "TS12_RIN_9.9_cDNA", "TS12_RIN_9.9_2_cDNA", "TS11_RIN_9.8_cDNA", "TS12_RIN_8.2_cDNA", "TS12_RIN_8.2_2_cDNA", "TS10_RIN_7.7_cDNA"
  )
  
  #remove rows that sum to zero
  CountMatrix <- CountMatrix[ rowSums(CountMatrix[,-1]) > 0, ]
  
  
  
  message("Importing gene length")
  lengths <- read.csv("/Users/yairp/Documents/Ph.D/Genomes/gfftools.gencode.v31_lengths.csv", header=T)
  #lengths <- read.csv("/Users/yairp/Documents/Ph.D/Genomes/Isoform.lengths_from_NananoCount.csv", header=T)
  
  txltmp <- as.data.frame(CountMatrix) %>% 
    tibble::rownames_to_column("Gene.ID") %>%
    tidyr::gather(key = sample, value = abundance, -Gene.ID) %>%
    tidyr::separate(sample, into = c("sample"), sep = "__")
    #dplyr::mutate(Gene.ID = gsub("\\.[^.]*$", "", Gene.ID)) # to remove the dots after the ENSG 
  
  # Add gene length here and merge data 
  Merged_lengths <- as.data.frame(merge(txltmp,lengths, by.x="Gene.ID", by.y="gene", all.x=FALSE)) %>% 
    tidyr::drop_na() %>% 
    dplyr::filter(median <= 5250)
  
  #Group lengths based breaks 
  #500bp breaks
  Merged_lengths$bin_Lengths <- cut(Merged_lengths$median, breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 
                                              4000, 4500, 5001), right = FALSE, dig.lab = 5) 
  
  #750bp breaks
  Merged_lengths$bin_Lengths <- cut(Merged_lengths$median, breaks = c(0, 750, 1500, 2250, 3000, 3750, 4500, 5250   
                                                                 ), right = FALSE, dig.lab = 5) 
  
  #1000bp breaks
  Merged_lengths$bin_Lengths <- cut(Merged_lengths$median, breaks = c(0, 1000, 2000, 3000, 4000, 5000), right = TRUE, dig.lab = 5) 
  
  Merged_lengths <- Merged_lengths %>% dplyr::mutate("sample" = 
                                                       case_when(sample == "TS12_RIN_9.9"  ~ 'CONTROL', 
                                                                 sample == "TS10_RIN_9.8" ~ 'CONTROL',
                                                                 sample == "TS11_RIN_9.7" ~ 'CONTROL',
                                                                 sample == "TS10_RIN_9.6" ~ 'VERY.HIGH',
                                                                 sample == "TS11_RIN_9.6" ~ 'VERY.HIGH',
                                                                 sample == "TS10_RIN_9.3" ~ 'HIGH',
                                                                 sample == "TS11_RIN_9.3" ~ 'HIGH',
                                                                 sample == "TS11_RIN_8.9" ~ 'MID',
                                                                 sample == "TS11_RIN_8.8" ~ 'MID',
                                                                 sample == "TS10_RIN_8.4" ~ 'MID',
                                                                 sample == "TS11_RIN_8.7" ~ 'LOW',
                                                                 sample == "TS12_RIN_8.2" ~ 'LOW',
                                                                 sample == "TS10_RIN_7.7" ~ 'LOW',
                                                                 sample == "TS12_RIN_7.3" ~ 'VERY.LOW',
                                                                 sample == "TS12_RIN_7.2" ~ 'VERY.LOW',
                                                                 sample == "TS12_RIN_9.9_cDNA"  ~ 'CONTROL.cDNA',
                                                                 sample == "TS12_RIN_9.9_2_cDNA"  ~ 'CONTROL.cDNA',
                                                                 sample == "TS11_RIN_9.8_cDNA"  ~ 'CONTROL.cDNA',
                                                                 sample == "TS12_RIN_8.2_cDNA"  ~ 'LOW.cDNA',
                                                                 sample == "TS12_RIN_8.2_2_cDNA"  ~ 'LOW.cDNA',
                                                                 sample == "TS10_RIN_7.7_cDNA"  ~ 'LOW.cDNA'))
  
  message("adding Time")
  Merged_lengths <- add_time(Merged_lengths) 
  
  message("adding Method")
  Merged_lengths <- add_method(Merged_lengths) 
  
  message("calculating detection rate")
  
  #Use only dRNA
  dRNA <- filter(Merged_lengths, Method  == 'dRNA')
  
  #need to turn this into a  function
  plot.me_lengths_1 <- dRNA %>% 
    group_by(Time, bin_Lengths) %>% 
    summarise(sum_var1 = sum(abundance)) %>%
    group_by(Time) %>%
    mutate(detRate = sum_var1/sum(sum_var1))
  
  #compare DRNA to cDNA
  dRNAvcDNA_strat <- filter(Merged_lengths, Method  == 'cDNA' | sample  == 'CONTROL' | sample  == 'LOW')
  plot.me_lengths_cdNA_DRNA <- dRNAvcDNA_strat %>% 
    group_by(Method, Time, bin_Lengths) %>% 
    summarise(sum_var1 = sum(abundance)) %>%
    group_by(Time, Method) %>%
    mutate(detRate = sum_var1/sum(sum_var1))
  
  
  message("Plotting Lengths vs detection rate") #proportion 
  
  pdf("Det_rate_cDNA_vsDRNA_Genes_750bp_barplot)by length.pdf", width=8, height=4) # run this for now but need to get rid for clean function
  ggplot(plot.me_lengths_cdNA_DRNA, aes(x = bin_Lengths , y = detRate)) + 
    geom_bar(aes(fill =Method),
             stat = "identity", position = position_dodge(0.85), width=0.8) + 
    theme_bw() + 
    #ylim(0,0.6) +
    #scale_fill_manual(values = ds_colors, name = "") + 
    facet_wrap(~ Time) + 
    xlab("Gene length") + ylab("Detected Gene proportion") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle("Genes Detected stratified by legnth (bp)") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
      plot.title = element_text(color="black", size=20),
      axis.title.x = element_text(color="black", size=14),
      axis.title.y = element_text(color="black", size=14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()) +
    theme(legend.title = element_blank())
  dev.off()
  
  
  pdf("Genes_stratified_length_750bp_barplot.pdf", width=8, height=4) # run this for now but need to get rid for clean function
  ggplot(plot.me_lengths_1, aes(x = bin_Lengths , y = detRate)) + 
    geom_bar(aes(fill =Time),
             stat = "identity", position = position_dodge(0.85), width=0.8) + 
    theme_bw() + 
    ylim(0,0.6) +
    #scale_fill_manual(values = ds_colors, name = "") + 
    #facet_wrap(~ Time) + 
    xlab("Gene length") + ylab("Detected Gene proportion") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle("Genes Detected stratified by legnth (bp)") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
      plot.title = element_text(color="black", size=20),
      axis.title.x = element_text(color="black", size=14),
      axis.title.y = element_text(color="black", size=14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()) +
    theme(legend.title = element_blank())
  dev.off()
    


}  


suppressWarnings(
  main())


setwd("~/Documents/Ph.D/RNA_degradation_study/Deseq_15_samples/20210729_ft_counts/29_July_feature_counts/")

#write.csv(plot.me_lengths_1, "1000bp_Detection_curves_forplotting_genes_median_length.csv" )

plot.me_lengths_1 <- read.csv("../1000bp_Detection_curves_forplotting_genes_median_length.csv")

#plot.me_lengths_1$bin_Lengths <- factor(plot.me_lengths_1$bin_Lengths, levels = c('[0,750)','[750,1500)','[1500,2250)','[2250,3000)','[3000,3750)','[3750,4500)', '[4500,5250)'))

plot.me_lengths_1$bin_Lengths <- factor(plot.me_lengths_1$bin_Lengths, levels = c('(0,1000]','(1000,2000]','(2000,3000]','(3000,4000]','(4000,5000]'))




#thought we could drawing a line graph instead. Might be a bit clearer 
pdf("1kb_Genes_Line_detection_curve_Normalised to control_28_4.pdf", width=8, height=4)
ggplot(plot.me_lengths_1,
       aes(x = bin_Lengths, y =detRate, colour=Time, group = Time)) +
  geom_line(size=0.5) + 
  geom_point(size=1) +
  geom_hline(yintercept=0, linetype="dashed", colour = "red") +
  theme_bw() +
  #ylim(-0.4,0.2) +
  #facet_wrap(~ Time) + 
  xlab("Gene length") + ylab("Detected Gene proportion") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Genes Detected stratified by legnth (bp)") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

ggplot(lengths) +
  geom_boxplot(aes(x=lengths$Length)) +
  xlim(0,50000)
summary(lengths)  



###########################################
#Stacked bar graph like 500bp bins test
#might be a little easier to visulaise 
##########################################

plot.me_lengths_1$bin_Lengths <- factor(plot.me_lengths_1$bin_Lengths, levels = c('[4500,5250)','[3750,4500)','[3000,3750)','[2250,3000)','[1500,2250)','[750,1500)','[0,750)'))

pdf("Genes_stratified_length_750bp_barplot_side_by_side.pdf", width=8, height=8) # run this for now but need to get rid for clean function
ggplot(plot.me_lengths_1, aes(x = detRate , y = Time)) + 
  geom_bar(aes(fill =bin_Lengths),
           stat = "identity", width = 0.7, position = 'fill') + 
  theme_bw() + 
  #scale_fill_manual(values = ds_colors, name = "") + 
  #facet_wrap(~ Time) + 
  xlab("Detected Gene proportion") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Genes Detected stratified by legnth (bp)") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    plot.title = element_text(color="black", size=20),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  theme(legend.title = element_blank())
dev.off()
