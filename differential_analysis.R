rm(list = ls())

#### 叶片中基因型1处理间比较
df1 <- read.table("G1T12L_counts.txt",header = T,stringsAsFactors = F,row.names = 1)
df1 <- as.matrix(df1)
glimpse(df1)
coldata <- data.frame(a = colnames(df1))
coldata$condition <- c(rep("treated1",3),rep("treated2",3))
rownames(coldata) <- coldata[,1]

coldata <- as.matrix(coldata)
coldata <- coldata[,-1]
coldata <- as.data.frame(coldata)
colnames(coldata) <- "condition"

dds <- DESeqDataSetFromMatrix(countData = df1,
                              colData = coldata,
                              design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition <- factor(dds$condition,levels = c("treated2","treated1"))

dds <- DESeq(dds)
res <- results(dds)
summary(res)
res05 <- results(dds,alpha = 0.05)
res051 <- as.data.frame(res05)
res051 %>% arrange(padj,log2FoldChange)

up <- subset(res051,log2FoldChange > 0 & padj < 0.05)


low <- subset(res051,log2FoldChange < 0 & padj < 0.05)
up <- res051 %>% filter(log2FoldChange > 0,padj < 0.05)
low <- res051 %>% filter(log2FoldChange < 0, padj < 0.05)
write.table(up,"G1_L_up.txt",sep = "\t",quote = F)
write.table(low,"G1_L_low.txt",sep = "\t",quote = F)

up_id <- rownames(up)
up_id <- data.frame(gene_id = up_id)

low_id <- rownames(low)
low_id <- data.frame(gene_id = low_id)


write.table(up_id,"G1_L_up_id.txt",quote = F)
write.table(low_id,"G1_L_low_id.txt",quote = F)
summary(res05)

plotMA(res01)

#热图
library(circlize)
ntd <- normTransform(dds)
mat <- assay(ntd)
low_res <- res051[rownames(res051) %in% low_id,]
low_res <- low_res[order(low_res$log2FoldChange,decreasing = F),]

low_resid <- rownames(low_res)

up_res <- res051[rownames(res051) %in% up_id,]
up_res <- up_res[order(up_res$log2FoldChange,decreasing = T),]
head(up_res)
up_resid <- rownames(up_res)
df5 <- read.table("G1T12L_counts.txt",header = T,row.names = 1)
df4 <- df5[rownames(df5) %in%  c(up_resid,low_resid),]
class(df4)
dim(df4)



dd <- t(apply(df4,1,scale))
colnames(dd) <- colnames(df4)
names(df4)

col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))
p1 <- Heatmap(dd,cluster_rows = T,cluster_columns = T,name = "",
        col = col_fun,
        column_title_side = "top",
        show_row_names = F
        )

df5 <- read.table("G1T12L_counts.txt",header = T,row.names = 1)


##plotMA 图
resultsNames(dds)
resLFC <- lfcShrink(dds,coef = "condition_treated1_vs_treated2")
plotMA(resLFC) 


#### 叶片中基因型2处理间比较

df2 <- read.table("G2T12L_counts.txt",header = T,stringsAsFactors = F,row.names = 1)
df2 <- as.matrix(df2)

glimpse(df2)
coldata <- data.frame(a = colnames(df2))
coldata$condition <- c(rep("treated1",3),rep("treated2",3))
rownames(coldata) <- coldata[,1]

coldata <- as.matrix(coldata)
coldata <- coldata[,-1]
coldata <- as.data.frame(coldata)
colnames(coldata) <- "condition"

dds <- DESeqDataSetFromMatrix(countData = df2,
                              colData = coldata,
                              design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition <- factor(dds$condition,levels = c("treated2","treated1"))

dds <- DESeq(dds)
res <- results(dds)
summary(res)
res05 <- results(dds,alpha = 0.05)
res051 <- as.data.frame(res05)
summary(res05)

up <- subset(res051,log2FoldChange > 0 & padj < 0.05)
low <- subset(res051,log2FoldChange < 0 & padj < 0.05)

write.table(up,"G2_L_up.txt",sep = "\t",quote = F)
write.table(low,"G2_L_low.txt",sep = "\t",quote = F)

up_id <- rownames(up)
up_id <- data.frame(gene_id = up_id)

low_id <- rownames(low)
low_id <- data.frame(gene_id = low_id)

write.table(up_id,"G2_L_up_id.txt",quote = F)
write.table(low_id,"G2_L_low_id.txt",quote = F)
summary(res05)

plotMA(res05)


#热图
library(circlize)

low_res <- res051[rownames(res051) %in% low_id,]
low_res <- low_res[order(low_res$log2FoldChange,decreasing = F),]
head(low_res)

low_resid <- rownames(low_res)


up_res <- up_res[order(up_res$log2FoldChange,decreasing = T),]
head(up_res)
up_resid <- rownames(up_res)

df5 <- read.table("G2T12L_counts.txt",row.names = 1,header = T)
df4 <- df5[rownames(df5) %in%  c(up_resid,low_resid),]
dim(df4)
head(df4)
class(df4)




dd <- t(apply(df4,1,scale))
colnames(dd) <- colnames(df4)


col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))
p2 <- Heatmap(dd,cluster_rows = T,cluster_columns = T,name = "",
              col = col_fun,
              column_title_side = "bottom",
              show_row_names = F
)

df5 <- read.table("G1T12L_counts.txt",header = T,row.names = 1)



##plotMA 图
resultsNames(dds)
resLFC <- lfcShrink(dds,coef = "condition_treated1_vs_treated2")
plotMA(resLFC)


##根系基因型1间的比较

df2 <- read.table("G1T12R_counts.txt",header = T,stringsAsFactors = F,row.names = 1)
df2 <- as.matrix(df2)

glimpse(df2)
coldata <- data.frame(a = colnames(df2))
coldata$condition <- c(rep("treated1",3),rep("treated2",3))
rownames(coldata) <- coldata[,1]

coldata <- as.matrix(coldata)
coldata <- coldata[,-1]
coldata <- as.data.frame(coldata)
colnames(coldata) <- "condition"

dds <- DESeqDataSetFromMatrix(countData = df2,
                              colData = coldata,
                              design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition <- factor(dds$condition,levels = c("treated2","treated1"))

dds <- DESeq(dds)
res <- results(dds)
summary(res)
res05 <- results(dds,alpha = 0.1)
res051 <- as.data.frame(res05)
summary(res05)

up <- subset(res05,log2FoldChange > 0 & padj < 0.1)
low <- subset(res05,log2FoldChange < 0 & padj < 0.1)

write.table(up,"G1_R_up.txt",sep = "\t",quote = F)
write.table(low,"G2_R_low.txt",sep = "\t",quote = F)

up_id <- rownames(up)
up_id <- data.frame(gene_id = up_id)

low_id <- rownames(low)
low_id <- data.frame(gene_id = low_id)

write.table(up_id,"G2_L_up_id.txt",quote = F)
write.table(low_id,"G2_L_low_id.txt",quote = F)
summary(res05)

# 热图
low_res <- res051[rownames(res051) %in% low_id,]
dim(low_res)
low_res <- low_res[order(low_res$log2FoldChange,decreasing = F),]
head(low_res)

low_resid <- rownames(low_res)

up_res <- res051[rownames(res051) %in% up_id,]
up_res <- up_res[order(up_res$log2FoldChange,decreasing = T),]
head(up_res)
up_resid <- rownames(up_res)

df5 <- read.table("G1T12R_counts.txt",row.names = 1,header = T)
df4 <- df5[rownames(df5) %in%  c(up_resid,low_resid),]
dim(df4)
head(df4)
class(df4)




dd <- t(apply(df4,1,scale))
colnames(dd) <- colnames(df4)


col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))
p2 <- Heatmap(dd,name = "",
              col = col_fun,
              column_title_side = "bottom",
              show_row_names = F
)

Heatmap
dim(dd)


##plotMA 图

#resultsNames(dds)
#resLFC <- lfcShrink(dds,coef = "condition_treated1_vs_treated2")
#plotMA(resLFC)

plotMA(res05)


##根系中基因型2处理间比较

df2 <- read.table("G2T12R_counts.txt",header = T,stringsAsFactors = F,row.names = 1)
df2 <- as.matrix(df2)

glimpse(df2)
coldata <- data.frame(a = colnames(df2))
coldata$condition <- c(rep("treated1",3),rep("treated2",3))
rownames(coldata) <- coldata[,1]

coldata <- as.matrix(coldata)
coldata <- coldata[,-1]
coldata <- as.data.frame(coldata)
colnames(coldata) <- "condition"

dds <- DESeqDataSetFromMatrix(countData = df2,
                              colData = coldata,
                              design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition <- factor(dds$condition,levels = c("treated2","treated1"))

dds <- DESeq(dds)
res <- results(dds)
summary(res)
res05 <- results(dds,alpha = 0.05)
res051 <- as.data.frame(res05)
summary(res05)

up <- subset(res051,log2FoldChange > 0 & padj < 0.1)
low <- subset(res051,log2FoldChange < 0 & padj < 0.1)

write.table(up,"G2_R_up.txt",sep = "\t",quote = F)
write.table(low,"G2_R_low.txt",sep = "\t",quote = F)

up_id <- rownames(up)
up_id <- data.frame(gene_id = up_id)

low_id <- rownames(low)
low_id <- data.frame(gene_id = low_id)

write.table(up_id,"G2_R_up_id.txt",quote = F)
write.table(low_id,"G2_R_low_id.txt",quote = F)
summary(res05)

# 热图
low_res <- res051[rownames(res051) %in% low_id,]
dim(low_res)
low_res <- low_res[order(low_res$log2FoldChange,decreasing = F),]
head(low_res)

low_resid <- rownames(low_res)

up_res <- res051[rownames(res051) %in% up_id,]
up_res <- up_res[order(up_res$log2FoldChange,decreasing = T),]
head(up_res)
up_resid <- rownames(up_res)

df5 <- read.table("G2T12R_counts.txt",row.names = 1,header = T)
df4 <- df5[rownames(df5) %in%  c(up_resid,low_resid),]
dim(df4)
head(df4)
class(df4)




dd <- t(apply(df4,1,scale))
colnames(dd) <- colnames(df4)


col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))
p2 <- Heatmap(dd,name = "",
              col = col_fun,
              column_title_side = "bottom",
              show_row_names = F
)





##plotMA 图

#resultsNames(dds)
#resLFC <- lfcShrink(dds,coef = "condition_treated1_vs_treated2")
#plotMA(resLFC)

plotMA(res05)

