setwd("/home/chenzhi/fareeha/count/results/clean_count/")

library(tidyverse)
library(DESeq2)

df <- read.table("all_counts.txt",header = T,sep = "\t",stringsAsFactors = F)

df %>% glimpse()
name <- names(df)[-1]
coldata <- data.frame(a = name,
                      condition = rep(c(rep("treated1",3),
                                        rep("treated2",3)),
                                      4))
coldata_L <- coldata[c(7:12,19:24),] 
coldata_L$type = c(rep("genetype1",6),rep("genetype2",6))

coldata_R <- coldata[-c(7:12,19:24),]
coldata_R$type = c(rep("genetype1",6),rep("genetype2",6))


#### input count matrix (Leaf)

dfG1L <- read.table("G1T12L_counts.txt",header = T,stringsAsFactors = F)
dfG2L <- read.table("G2T12L_counts.txt",header = T,stringsAsFactors = F)

dfG12L <- left_join(dfG1L,dfG2L)
glimpse(dfG12L)
rownames(dfG12L) <- dfG12L[,1]
dfG12L <- dfG12L[,-1]
head(dfG12L)

dfG12L <- as.matrix(dfG12L)

rownames(coldata_L) <- coldata_L[,1]
coldata_L <- coldata_L[,-1]

dds <- DESeqDataSetFromMatrix(countData = dfG12L,
                              colData = coldata_L,
                              design = ~ condition)
dds

vsd <- vst(dds,blind = F)
head(assay(vsd),3)
plotPCA(vsd,intgroup = c("condition","type"))

pcaData <- plotPCA(vsd,intgroup = c("condition","type"),
                   returnData = T)
pcaData
percentVar <- round(100 * attr(pcaData,"percentVar"))

# make PCAplot
ggplot(pcaData,aes(PC1,PC2,color = condition,shape = type)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance" )) +
  ylab(paste0("PC2: ",percentVar[2],"% variance" )) +
  coord_fixed() + 
  theme_bw() +
  ggtitle("PCA of Leaf tissues") +
  theme(plot.title = element_text(hjust = 0.5))
  

#input count matrix (root)

dfG1R <- read.table("G1T12R_counts.txt",header = T,stringsAsFactors = F)
dfG2R <- read.table("G2T12R_counts.txt",header = T,stringsAsFactors = F)

dfG12R <- left_join(dfG1R,dfG2R)
dfG12R1 <- dfG12R[sample(nrow(dfG12R),15000),]
glimpse(dfG12R)
rownames(dfG12R1) <- dfG12R1[,1]
dfG12R1 <- dfG12R1[,-1]
head(dfG12R1)
dfG12R <- as.matrix(dfG12R1)

rownames(coldata_R) <- coldata_R[,1]
coldata_R <- coldata_R[,-1]

dds1 <- DESeqDataSetFromMatrix(countData = dfG12R,
                              colData = coldata_R,
                              design = ~ condition)
dds1

vsd1 <- vst(dds1,blind = F)
head(assay(vsd1),3)
plotPCA(vsd1,intgroup = c("condition","type"))

pcaData <- plotPCA(vsd1,intgroup = c("condition","type"),
                   returnData = T)
pcaData
percentVar <- round(100 * attr(pcaData,"percentVar"))

# make PCAplot
ggplot(pcaData,aes(PC1,PC2,color = condition,shape = type)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance" )) +
  ylab(paste0("PC2: ",percentVar[2],"% variance" )) +
  coord_fixed() + 
  theme_bw() +
  ggtitle("PCA of Root tissues") +
  theme(plot.title = element_text(hjust = 0.5))
