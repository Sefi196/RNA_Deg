library(DESeq2)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(ggplot2)
library(d3heatmap)
library(data.table)

#Can be run for either 5y or Sequin counts 
#Gene counts -> ft.counts
#Isoform counts -> NanoCount

#set the working dir
#For Genens
setwd("~/Documents/Ph.D/RNA_degradation_study/Deseq_15_samples/20210729_ft_counts/29_July_feature_counts/")
countdata <- read.csv("ft.counts.csv", row.names = 1, header = T)

#for Isoforms 
setwd("~/Documents/Ph.D/RNA_degradation_study/Deseq_15_samples/NanoCount-0.4.0/")
countdata <- read.csv("NanoCount.Count.Matrix.csv", row.names = 1, header = T)

#to make as matric instead of df
countdata <- as.matrix(countdata)
head(countdata)

#Make correct names
colnames(countdata) <- c("TS12_RIN_9.9", "TS10_RIN_9.8", "TS11_RIN_9.7", "TS10_RIN_9.6","TS11_RIN_9.6", "TS10_RIN_9.3", "TS11_RIN_9.3", "TS11_RIN_8.9", "TS11_RIN_8.8", "TS10_RIN_8.4", "TS11_RIN_8.7", "TS12_RIN_8.2", "TS10_RIN_7.7", "TS12_RIN_7.3", "TS12_RIN_7.2")

# Assign conditions
# This step is easy like below if samples from each condition are next to eachother in the matrix

# Normal condition 
Time <- as.factor(c(rep("0", 3), rep("0.5", 2), rep("1", 2), rep("3_4", 3), rep("6", 3),rep("8", 2)))
seq_batch <- as.factor(c(2, 1, 4, 1, 4, 1, 2, 4, 3, 3, 1, 3, 1, 2, 4))
RIN <- scale(c(9.9, 9.8, 9.7, 9.6, 9.6, 9.3, 9.3, 8.9, 8.8, 8.4, 8.7, 8.2, 7.7, 7.3, 7.2))
extraction_batch <- as.factor(c(3, 1, 2, 1 , 2, 1, 2, 2, 2, 1, 2, 3, 1, 3, 3))

#Sequin conditions 
#Time <- as.factor(c(rep("T1", 3), rep("T6", 5)))
#seq_batch <- as.factor(c(2, 1, 3, 1, 1, 2, 1, 1))
#RIN <- scale(c(9.9, 9.8, 9.7, 9.6, 9.3, 9.3, 8.4, 7.7))
#extraction_batch <- as.factor(c(3, 1, 2, 1, 1, 2, 1, 1))
#deseq1<-data.frame(Time,RIN,seq_batch,extraction_batch)

deseq1<-data.frame(Time,RIN,seq_batch,extraction_batch)

# Make DESeq dataset
(coldata <- data.frame(row.names=colnames(countdata), deseq1))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~seq_batch + extraction_batch + Time)

# Optional filtering step to remove very low counts
#keep <- rowSums(counts(dds)) >= 50
#dds <- dds[keep,]

#filtering step to step to adress convergfance problem (found in bioconductor support) # only used for transcripts 
# probaly too stringet as reads here evn at the gene count are low. 
#filtering here means nc has to be greater than or equal to 10 in atleast 4 of the sa mples
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 5
dds <- dds[filter,]

# running deseq with adjusted iterations to adress converagce issue 
#dds <- estimateSizeFactors(dds)
#dds <- estimateDispersions(dds)
#dds <- nbinomWaldTest(dds, maxit=500)

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

#now filter out the significant genes padj < 0.05

table(res_T_0h_vs_0.5$padj<0.05)
table(res_T_0h_vs_1$padj<0.05)
table(res_T_0h_vs_3_4$padj<0.05)
table(res_T_0h_vs_6$padj<0.05)
table(res_T_0h_vs_8$padj<0.05)

# Normal conditions
res <- DESeq2::results(dds)
res

res <- res[order(res$padj), ]

table(res$padj<0.05)

#try plotting counts

topGene <- rownames(res)[which.min(res$padj)] # without RIN in the model, the top gene is ENSG00000166598.15
plotCounts(dds, gene = "ENSG00000166598.15", intgroup=c("Time"), transform = TRUE)

for (i in 1:5) { 
geneCounts <- plotCounts(dds, gene = rownames(res)[i], intgroup = c("Time"),
                         returnData = TRUE, transform = TRUE)
print(ggplot(geneCounts, aes(x = Time, y = count, group = "Time")) + 
  geom_point() + geom_line() +
  stat_summary(fun=mean, geom="line") +
  scale_y_log10() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = i, size=5,
       colour = "Time (h)"))
}

ggplot(geneCounts, aes(x = Time, y = count, group = "Time")) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = i, size=5,
       colour = "Time (h)")

geneCounts <- plotCounts(dds, gene = rownames(res)[1], intgroup = c("Time"),
                         returnData = TRUE)

ggplot(geneCounts, aes(x = Time, y = count, group = "Time")) + 
        geom_point() + stat_summary(fun=mean, geom="line") +
        scale_y_log10()


# Save normalised counts separately 
norm_counts <- as.data.frame(counts(dds, normalized=TRUE))
boxplot(norm_counts)

#write.csv(norm_counts, file="norm_counts.csv")
table(res_T_0h_vs_0.5$padj<0.05)
table(res_T_0h_vs_1$padj<0.05)
table(res_T_0h_vs_3_4$padj<0.05)
table(res_T_0h_vs_6$padj<0.05)
table(res_T_0h_vs_8$padj<0.05)

# Order by adjusted p-value
res1_T <- res[order(res_T_0h_vs_0.5$padj), ]
res2_T <- res[order(res_T_0h_vs_1$padj), ]
res3_T <- res[order(res_T_0h_vs_3_4$padj), ]
res4_T <- res[order(res_T_0h_vs_6$padj), ]
res5_T <- res_T_0h_vs_8[order(res_T_0h_vs_8$log2FoldChange), ]

# Merge with normalized count data
resdata_TTR_0.5 <- merge(as.data.frame(res_T_0h_vs_0.5), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resdata_TTR_1 <- merge(as.data.frame(res_T_0h_vs_1), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resdata_TTR_3_4 <- merge(as.data.frame(res_T_0h_vs_3_4), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resdata_TTR_6 <- merge(as.data.frame(res_T_0h_vs_6), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resdata_TTR_8 <- merge(as.data.frame(res_T_0h_vs_8), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

#merge witb Norm counts Sequin
resdata_seq <- merge(as.data.frame(res_seq), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata_seq)[1] <- "Gene_ID"


names(resdata_TTR_0.5)[1] <- "Gene_ID"
names(resdata_TTR_1)[1] <- "Gene_ID"
names(resdata_TTR_3_4)[1] <- "Gene_ID"
names(resdata_TTR_6)[1] <- "Gene_ID"
names(resdata_TTR_8)[1] <- "Gene_ID"


head(resdata_TG_0.5)

summary(res1_T)

# Write DE results 
write.csv(resdata_TT_0.5, file="diff_transctipt_0h_vs_0.5h_Time.csv")
write.csv(resdata_TT_1, file="diff_transctipt_0h_vs_1h_Time.csv")
write.csv(resdata_TT_3_4, file="diff_transctipt_0h_vs_3-4h_Time.csv")
write.csv(resdata_TT_6, file="diff_transctipt_0h_vs_6h_Time.csv")
write.csv(resdata_TTR_8, file="diff_transctipt_0h_vs_8h_full_model_No_RIN.csv")

# Plot dispersions
plotDispEsts(dds, main="Dispersion plot", genecol = "black", fitcol = "cyan", finalcol = "blue", legend = TRUE)

# RLD for viewing
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Plot residual p-values
hist(res$pvalue, breaks=50, col="grey")

#Set colours for plotting
# for full set of all 15 samples
mycols <- brewer.pal(8, "Accent")[1:length(unique(Time))]

# Heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
pdf("heatmap-samples.pdf", 18, 18, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[seq_batch], RowSideColors=mycols[extraction_batch],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# PCA
rld_pca <- function (rld, intgroup = "Time", ntop = 500, colors=NULL, main="Principal Component Analysis", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1) # change the [X,this] to get other PC components 
  pc1lab <- paste0("PC1: ",as.character(pc1var),"% variance")
  pc2lab <- paste0("PC2: ",as.character(pc2var),"% variance")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx, offset=1)) #makes sample labells
  #legend(legendpos, legend=levels(fac), col=colors, pch=6)
  legend("topright", legend=levels(fac), col=colors, pch = 6,
         xpd=TRUE, horiz=FALSE, bty="n"
  )
}
pdf("pca-genes 5y_1-2.pdf", 6, 6, pointsize=13)
rld_pca(rld, colors=mycols, intgroup="Time", xlim=c(-20, 20), ylim=c(-10, 10))
dev.off()

# MA Plot
maplot <- function (res, thresh=0.05, labelsig=FALSE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="blue", pch=20, cex=1))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
pdf("diffexpr-maplot-0.05.pdf", 18, 18, pointsize=20)
maplot(resdata_TTR_8, main="MA Plot")
dev.off()

# Volcano Plot
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, xlab="log2(Fold Change)", legendpos="topright", labelsig=FALSE, textcx=1.5, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, xlab=xlab, cex.axis=1.8, cex.lab=1.5, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("p-adj<",sigthresh,sep=""), paste("|log2(FC)|>",lfcthresh,sep=""), "both"), cex=1.5, pch=20, col=c("blue","orange","green"))
}
pdf("diffexpr-volcanoplot-hi-res.pdf", 18, 18, pointsize=20)
volcanoplot(resdata, lfcthresh=2, sigthresh=0.05, xlim=c(-7, 7), ylim=c(0,20), legendpos="topright")
dev.off()
