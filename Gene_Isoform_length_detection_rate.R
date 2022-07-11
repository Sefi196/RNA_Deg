library(ggplot2)
library(dplyr)

#Count matrix can be Genes or isofrom and/or cDNA. 

#set working directory 
setwd("")

#input count matrix
CountMatrix <- read.csv("~/Documents/Ph.D/RNA_degradation_study/Deseq_15_samples/20210729_ft_counts/29_July_feature_counts/ft.counts.csv", row.names = 1, header = T)

#CountMatrix <- read.csv("PCS110_dRNA_ft.counts.csv", row.names = 1, header = T)

CountMatrix <- as.matrix(as.data.frame(CountMatrix))
colnames(CountMatrix) <- c("TS12_RIN_9.9", "TS10_RIN_9.8", "TS11_RIN_9.7", "TS10_RIN_9.6","TS11_RIN_9.6", "TS10_RIN_9.3", "TS11_RIN_9.3", "TS11_RIN_8.9", "TS11_RIN_8.8", 
                           "TS10_RIN_8.4", "TS11_RIN_8.7", "TS12_RIN_8.2", "TS10_RIN_7.7", "TS12_RIN_7.3", "TS12_RIN_7.2")
                           #, "TS12_RIN_9.9_cDNA", "TS12_RIN_9.9_2_cDNA", "TS11_RIN_9.8_cDNA", "TS12_RIN_8.2_cDNA", "TS12_RIN_8.2_2_cDNA", "TS10_RIN_7.7_cDNA")

#remove rows that sum to zero
CountMatrix <- CountMatrix[ rowSums(CountMatrix[,-1]) > 0, ]


message("Importing gene length")
lengths <- read.csv("/lengths_GTFtools.gencode.v31.basic.annotation.csv", header=T)

txltmp <- as.data.frame(CountMatrix) %>% 
  tibble::rownames_to_column("Gene.ID") %>%
  tidyr::gather(key = sample, value = abundance, -Gene.ID) %>%
  tidyr::separate(sample, into = c("sample"), sep = "__")
#dplyr::mutate(Gene.ID = gsub("\\.[^.]*$", "", Gene.ID)) # to remove the dots after the ENSG 

# Add gene length here and merge data 
Merged_lengths <- as.data.frame(merge(txltmp,lengths, by.x="Gene.ID", by.y="gene", all.x=FALSE)) %>% 
  tidyr::drop_na() #%>% 
  #dplyr::filter(median <= 5000)

#Group lengths based breaks 
#500bp breaks
#Merged_lengths$bin_Lengths <- cut(Merged_lengths$median, breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 
#                                                                   4000, 4500, 5001), right = FALSE, dig.lab = 5) 

#750bp breaks
#Merged_lengths$bin_Lengths <- cut(Merged_lengths$median, breaks = c(0, 750, 1500, 2250, 3000, 3750, 4500, 5250   
#), right = FALSE, dig.lab = 5) 

#1000bp breaks
Merged_lengths$bin_Lengths <- cut(Merged_lengths$median, breaks = c(0, 1000, 2000, 3000, 4000, 500000), right = TRUE, dig.lab = 5) 

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

message("adding Time")
Merged_lengths <- add_time(Merged_lengths) 

message("adding Method")
Merged_lengths <- add_method(Merged_lengths) 

message("calculating detection rate")

#Use only dRNA
dRNA <- filter(Merged_lengths, Method  == 'dRNA')

#need to turn this into a  function
plot.me_lengths_genbasic <- dRNA %>% 
  group_by(Time, bin_Lengths) %>% 
  summarise(sum_var1 = sum(abundance)) %>%
  group_by(Time) %>%
  mutate(detRate = sum_var1/sum(sum_var1))

#Run if comapring DRNA to cDNA 
dRNAvcDNA_strat <- filter(Merged_lengths, Method  == 'cDNA' | sample  == 'CONTROL' | sample  == 'LOW')
plot.me_lengths_cdNA_DRNA <- dRNAvcDNA_strat %>% 
  group_by(Method, Time, bin_Lengths) %>% 
  summarise(sum_var1 = sum(abundance)) %>%
  group_by(Time, Method) %>%
  mutate(detRate = sum_var1/sum(sum_var1))


#####
Plotting Lengths vs detection rate
#######


pdf("Genes_stratified_length_750bp_barplot.pdf", width=8, height=4)
ggplot(plot.me_lengths_genbasic, aes(x = bin_Lengths , y = detRate)) + 
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

#Plot line graph 
pdf("1kb_Genes_Line_detection_curve_basic_gencode_gene_set_27_5_22.pdf", width=8, height=4)
ggplot(plot.me_lengths_genbasic,
       aes(x = bin_Lengths, y =detRate, colour=Time, group = Time)) +
  geom_line(size=0.5) + 
  geom_point(size=1) +
  geom_hline(yintercept=0, linetype="dashed", colour = "red") +
  theme_bw() +
  ylim(0,0.6) +
  #facet_wrap(~ Time) + 
  xlab("Gene length") + ylab("Detected Gene proportion") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Genes Detected stratified by legnth (bp)") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#write data frame to working dir
write.csv(plot.me_lengths_genbasic, "Basic_1000bp_Detection_curves_forplotting_genes_median_length.csv" )

####
#Run if comapring DRNA to cDNA 
####

pdf("Det_rate_cDNA_vsDRNA_Genes_750bp_barplot)by length.pdf", width=8, height=4) 
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



