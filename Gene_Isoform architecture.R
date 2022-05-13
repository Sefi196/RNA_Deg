library(data.table)
library(DESeq2)
library(dplyr)
library(reshape2)
library(tibble)
library(tidyr)

############
#GENE / Isoform architecture  
############

#Inputs 
# Count Matrix (Gene or Isoform)
# gtf File 
# biomart info of interest - e.g GC content, etc.

# Import Count data (can be from ft counts or NanoCount)  
countdata <- read.csv("ft.counts.csv", row.names = 1, header = T)

#Label col names 
colnames(countdata) <- c("TS12_RIN_9.9", "TS10_RIN_9.8", "TS11_RIN_9.7", "TS10_RIN_9.6","TS11_RIN_9.6", "TS10_RIN_9.3", "TS11_RIN_9.3", "TS11_RIN_8.9", "TS11_RIN_8.8", "TS10_RIN_8.4", "TS11_RIN_8.7", "TS12_RIN_8.2", "TS10_RIN_7.7", "TS12_RIN_7.3", "TS12_RIN_7.2")

#Set up de-seq pipeline to normalise counts 
# Assign conditions
Time <- as.factor(c(rep("0", 3), rep("0.5", 2), rep("1", 2), rep("3_4", 3), rep("6", 3),rep("8", 2)))
seq_batch <- as.factor(c(2, 1, 4, 1, 4, 1, 2, 4, 3, 3, 1, 3, 1, 2, 4))
RIN <- scale(c(9.9, 9.8, 9.7, 9.6, 9.6, 9.3, 9.3, 8.9, 8.8, 8.4, 8.7, 8.2, 7.7, 7.3, 7.2))
extraction_batch <- as.factor(c(3, 1, 2, 1 , 2, 1, 2, 2, 2, 1, 2, 3, 1, 3, 3))

deseq1<-data.frame(Time,RIN,seq_batch,extraction_batch)

# Make DESeq dataset
(coldata <- data.frame(row.names=colnames(countdata), deseq1))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~seq_batch + extraction_batch + Time)

#filtering here means nc has to be greater than or equal to 10 in at least 5 of the samples
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 5
dds <- dds[filter,]

# Run DESeq2 pipeline
dds <- DESeq(dds)

#check the levels 
resultsNames(dds)

# will try and use contrast condition
res_T_0h_vs_0.5 <- as.data.frame(results(dds, contrast = c("Time","0.5","0")))
res_T_0h_vs_1 <- as.data.frame(results(dds, contrast = c("Time","1","0")))
res_T_0h_vs_3_4 <- as.data.frame(results(dds, contrast = c("Time","3_4","0")))
res_T_0h_vs_6 <- as.data.frame(results(dds, contrast = c("Time","6","0")))
res_T_0h_vs_8 <- as.data.frame(results(dds, contrast = c("Time","8","0")))


############
#Cluster based on logFC 
############

res_T_0h_vs_0.5$Tm <- "0.5"
res_T_0h_vs_1$Tm <- "1"
res_T_0h_vs_3_4$Tm <- "3"
res_T_0h_vs_6$Tm <- "6"
res_T_0h_vs_8$Tm <- "8"

#Define function for renaming cols
Gene_names <- function(df) {
  df <- as.data.frame(df)
  df <- setDT(df, keep.rownames = TRUE)[]
  names(df)[1] <- "Gene.ID"
  return(df)
}

res_T_0h_vs_0.5 <- Gene_names(res_T_0h_vs_0.5)
res_T_0h_vs_1 <- Gene_names(res_T_0h_vs_1)
res_T_0h_vs_3_4 <- Gene_names(res_T_0h_vs_3_4)
res_T_0h_vs_6 <- Gene_names(res_T_0h_vs_6)
res_T_0h_vs_8 <- Gene_names(res_T_0h_vs_8)


names(res_T_0h_vs_0.5)[3] <- "logFC_0.5"
names(res_T_0h_vs_1)[3] <- "logFC_1"
names(res_T_0h_vs_3_4)[3] <- "logFC_3"
names(res_T_0h_vs_6)[3] <- "logFC_6"
names(res_T_0h_vs_8)[3] <- "logFC_8"

#Define function for extracting cols of interest 
pull_col <- function(df) {
  df[, c(1, 3)]
}

x <- pull_col(res_T_0h_vs_0.5)
y <- pull_col(res_T_0h_vs_1)
z <- pull_col(res_T_0h_vs_3_4)
xx <- pull_col(res_T_0h_vs_6)
yy <- pull_col(res_T_0h_vs_8)

#Combine data frames together and add LF change 0 at time=0
LFC <- Reduce(inner_join, list(x, y, z, xx, yy))
LFC <- LFC %>% remove_rownames %>% column_to_rownames(var="Gene.ID")
LFC$logFC_0 <- 0 

#Move the last col to the first col
LFC <- LFC %>%
  dplyr::select(logFC_0, everything())

#run K means clustering function on LFC data frame.  
Kmean_LFC <- kmeans(LFC, 4, iter.max = 100, nstart = 10)
Kmean_LFC$size
kClusters <- Kmean_LFC$cluster
kCluster.df <- Gene_names(kClusters) 


#Now we can calculate the cluster ‘cores’ aka centroids:
# function to find centroid in cluster i
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}

kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, LFC, kClusters)

#Plotting the centroids to see how they behave:
#get in long form for plotting
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')

p1 <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Time",color = "Cluster")
p1
