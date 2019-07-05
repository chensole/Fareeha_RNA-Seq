library(clusterProfiler)

# G1T12L GO 
DElist_low <- read.table("../G1_L/G1_L_low_id.txt")
DElist_up <- read.table("../G1_L/G1_L_up_id.txt")
colnames(DElist_up) <- "gene_id"
id <- rbind(DElist_low,DElist_up)
id$gene_id
term2gene <- read.table("../../../../ref/term_gene.txt",header = F)
term2name <- read.csv("../../../../ref/term_name.txt",sep = "\t",header = F)
dim(term2name)
table(unique(term2gene$V1) %in% term2name$V1)

x <- enricher(id$gene_id,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = 0.01,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.05)
barplot(x)

out <- as.data.frame(x)

write.table(out,"G1T12L_GO_enrichment_results.txt",sep = "\t",quote = F)
getwd()



# G1T12R
DElist_low <- read.table("../G1_R/G1_R_low_id.txt")
DElist_up <- read.table("../G1_R/G1_R_up_id.txt")
colnames(DElist_up) <- "gene_id"
id <- rbind(DElist_low,DElist_up)
id$gene_id
term2gene <- read.table("../../../../ref/term_gene.txt",header = F)
term2name <- read.csv("../../../../ref/term_name.txt",sep = "\t",header = F)
dim(term2name)
table(unique(term2gene$V1) %in% term2name$V1)

x <- enricher(id$gene_id,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = 0.01,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.05)
barplot(x)
dotplot(x)
out <- as.data.frame(x)

write.table(out,"G1T12R_GO_enrichment_results.txt",sep = "\t",quote = F)


# G2T12L GO

DElist_low <- read.table("../G2_L/G2_L_low_id.txt")
DElist_up <- read.table("../G2_L/G2_L_up_id.txt")
colnames(DElist_up) <- "gene_id"
id <- rbind(DElist_low,DElist_up)
id$gene_id
term2gene <- read.table("../../../../ref/term_gene.txt",header = F)
term2name <- read.csv("../../../../ref/term_name.txt",sep = "\t",header = F)
dim(term2name)
table(unique(term2gene$V1) %in% term2name$V1)

x <- enricher(id$gene_id,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = 0.01,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.05)
barplot(x)
dotplot(x)
out <- as.data.frame(x)

write.table(out,"G2T12L_GO_enrichment_results.txt",sep = "\t",quote = F)


# G2T12R GO

# G2T12L GO

DElist_low <- read.table("../G2_R/G2_R_low_id.txt")
DElist_up <- read.table("../G2_R/G2_R_up_id.txt")
colnames(DElist_up) <- "gene_id"
id <- rbind(DElist_low,DElist_up)
id$gene_id
term2gene <- read.table("../../../../ref/term_gene.txt",header = F)
term2name <- read.csv("../../../../ref/term_name.txt",sep = "\t",header = F)
dim(term2name)
table(unique(term2gene$V1) %in% term2name$V1)

x <- enricher(id$gene_id,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = 0.01,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.05)
barplot(x)
dotplot(x)
out <- as.data.frame(x)

write.table(out,"G2T12R_GO_enrichment_results.txt",sep = "\t",quote = F)
