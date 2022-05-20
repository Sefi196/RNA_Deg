#Inputs are count Matrix for genes or isoforms

library(DESeq2)
library(data.table)
library(dplyr)
library(ggplot2)

# Import data 
countdata <- read.csv("ft.counts.csv", row.names = 1, header = T)

#to make as matrix instead of df
countdata <- as.data.frame(countdata)
head(countdata)

#Make correct names
colnames(countdata) <- c("TS12_RIN_9.9", "TS10_RIN_9.8", "TS11_RIN_9.7", "TS10_RIN_9.6","TS11_RIN_9.6", "TS10_RIN_9.3", "TS11_RIN_9.3", "TS11_RIN_8.9", "TS11_RIN_8.8", "TS10_RIN_8.4", "TS11_RIN_8.7", "TS12_RIN_8.2", "TS10_RIN_7.7", "TS12_RIN_7.3", "TS12_RIN_7.2")

# Normal condition 
Time <- as.factor(c(rep("0", 3), rep("0.5", 2), rep("1", 2), rep("3_4", 3), rep("6", 3),rep("8", 2)))
seq_batch <- as.factor(c(2, 1, 4, 1, 4, 1, 2, 4, 3, 3, 1, 3, 1, 2, 4))
RIN <- scale(c(9.9, 9.8, 9.7, 9.6, 9.6, 9.3, 9.3, 8.9, 8.8, 8.4, 8.7, 8.2, 7.7, 7.3, 7.2))
extraction_batch <- as.factor(c(3, 1, 2, 1 , 2, 1, 2, 2, 2, 1, 2, 3, 1, 3, 3))

deseq1<-data.frame(Time,RIN,seq_batch,extraction_batch)

#Use deseq to norm counts
coldata <- data.frame(row.names=colnames(countdata), deseq1)
dds_n <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~Time)
dds_n <- estimateSizeFactors(dds_n)
dds_n <- DESeq(dds_n)

Total_count <- as.data.frame(counts(dds_n, normalized=TRUE))

# can use nc counts from Here
# Remove rows with all zero values
not_zero <- Total_count %>% filter(across(everything(.)) != 0)

Pull out rows from NC counts that
#5 NC in at least 4 of the contro and 0.5 hours samples 
df.1 <- subset(not_zero, rowSums(not_zero[1:5] > 5) > 3) 

#0 NC in both 8 hour samples 
df.2 <- subset(df.1, rowSums(df.1[14:15] < 0.000001) > 1)

##############
plot these genes 
##############

#convert to plotting format 
df.3 <- setDT(df.2, keep.rownames = TRUE)[]
names(df.3)[1] <- "Gene.ID"
df.3molten <- melt(df.3)
colnames(df.3molten) <- c('Gene.ID','Sample','Norm_count')


# Add R2 values to a fitted model. 
#build data frame in the correct format. Make sure the first col is labelled Gene.ID
D <- as.data.frame(df.3)

#Converge data structure frame based on Gene ID 
D <- D %>% 
  gather(Sample, Count, -Gene.ID)
head(D)

#add log+1 cpount col -> will make plotting much easier 
D$log_Norm_count <- log(D$Count+1) 

#use mutate to add a new coll based on the the Sample ID - For me im adding Time in minutes. This will be the dependant variable  
D <- D %>% mutate(Time =
                    case_when(Sample == "TS12_RIN_9.9"  ~ 0, 
                              Sample == "TS10_RIN_9.8" ~ 0,
                              Sample == "TS11_RIN_9.7" ~ 0,
                              Sample == "TS10_RIN_9.6" ~ 0.5,
                              Sample == "TS11_RIN_9.6" ~ 0.5,
                              Sample == "TS10_RIN_9.3" ~ 1,
                              Sample == "TS11_RIN_9.3" ~ 1,
                              Sample == "TS11_RIN_8.9" ~ 3,
                              Sample == "TS11_RIN_8.8" ~ 4,
                              Sample == "TS10_RIN_8.4" ~ 3,
                              Sample == "TS11_RIN_8.7" ~ 6,
                              Sample == "TS12_RIN_8.2" ~ 6,
                              Sample == "TS10_RIN_7.7" ~ 6,
                              Sample == "TS12_RIN_7.3" ~ 8,
                              Sample == "TS12_RIN_7.2" ~ 8)
)

#following the worked example in R for data science. 
by_gene <- D %>% 
  group_by(Gene.ID) %>% 
  nest()

#to see the data that has been nested at gene 15 > also important to check that the col are of the correct format
by_gene$data[[10]]

#model fitting function - exponential decay
model_exp <- function(df) {
  lm(log(Count+1) ~ Time, data = df)
}

#model fitting function - linear Model
model_lin <- function(df) {
  lm(log_Norm_count ~ Time, data = df)
}

Emodel <- map(by_gene$data, model_exp)

by_gene <- by_gene %>% 
  mutate(model = map(data, model_exp))
by_gene

#Model Quality
#this adds R2 values to the table
by_gene <- by_gene %>% 
  mutate(glance = map(model, broom::glance)) %>% 
  unnest(glance, .drop = TRUE)

by_gene$r = sqrt(by_gene$r.squared)

#remove the poor fits from the plotme group
plotme_new <- filter(by_gene, r >0.5)

plotme_new <- plotme_new %>% unnest(data)

pdf("5counts_inat_least_5_samples_Genes_that_can't_be detected_after_8hours_test.pdf", width = 4, height=2) 
ggplot(plotme_new) +
  geom_point(aes(x=Time, y=Count), size =0.8) +
  geom_smooth(aes(x=Time, y=Count), method = "glm", method.args = list(family = poisson), colour = "RED", se = FALSE, size =0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)+
  facet_wrap(facets = 'Gene.ID') +
  ylab("Normalised Counts") +
  xlab("Time (min)") +
  #ylim(0,20) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    plot.title = element_text(color="black", size=20),
    axis.title.x = element_text(color="black", size=10),
    axis.title.y = element_text(color="black", size=10),
    panel.grid.major = element_blank(),
    strip.text = element_text(size=7))
dev.off()

pdf("5counts_inat_least_5_samples_Genes that dissapear_logCountsplus linear model.pdf", width = 8, height=8) 
ggplot(plotme_new) +
  geom_point(aes(x=Time, y=log_Norm_count), size =0.8) +
  #geom_smooth(aes(x=Time, y=Count), method = "glm", method.args = list(family = poisson), colour = "RED", se = FALSE, size =0.5) +
  geom_smooth(aes(x=Time, y=log_Norm_count), method = "lm", formula = y ~ x, colour = "RED", se = FALSE, size =0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5)+
  #facet_wrap(facets = 'Gene.ID', scales = "free") +
  facet_wrap(facets = 'Gene.ID') +
  #scale_y_continuous(breaks = integer_breaks()) +
  ylab("log(Normalised Counts + 1)") +
  xlab("Time (hours)") +
  #ylim(-5,70) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    plot.title = element_text(color="black", size=20),
    axis.title.x = element_text(color="black", size=10),
    axis.title.y = element_text(color="black", size=10),
    panel.grid.major = element_blank(),
    strip.text = element_text(size=7))
#panel.grid.minor = element_blank())
dev.off()

#remove nested model to save csv
df <- plotme_new[ -c(6) ]

write.csv(df, file="genes_dissapear_8hours_5counts_inat_least_5_samples.csv")
