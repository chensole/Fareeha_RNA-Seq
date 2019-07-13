library(clusterProfiler)
setwd("/home/chenzhi/fareeha/Results_to_fareeha")
pathway_gene <- read.table("../fareeha/ref/pathway_2_gene.txt",header = F)
pathway_name <- read.csv("../fareeha/ref/pathway_2_name.txt",header = F,sep = ";")



# G1T12L kegg
DElist_low <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype1_leaf_treat1_VS_treat2/G1_L_low_expression_id.txt")
DElist_up <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype1_leaf_treat1_VS_treat2/G1_L_up_expression_id.txt")
colnames(DElist_up) <- "gene_id"
id <- rbind(DElist_low,DElist_up)
id$gene_id



x <- enricher(id$gene_id,
              TERM2GENE = pathway_gene,
              TERM2NAME = pathway_name,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.1)
dotplot(x)

out <- as.data.frame(x)

write.table(out,"G1T12L_kegg_enrichment_results.txt",sep = "\t",quote = F)
getwd()



# G1T12R kegg
DElist_low <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype1_root_treat1_VS_treat2/G1_R_low_expression_id.txt")
DElist_up <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype1_root_treat1_VS_treat2/G1_R_up_expression_id.txt")
colnames(DElist_up) <- "gene_id"
id <- rbind(DElist_low,DElist_up)
id$gene_id



x <- enricher(id$gene_id,
              TERM2GENE = pathway_gene,
              TERM2NAME = pathway_name,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.1)

dotplot(x)
out <- as.data.frame(x)

write.table(out,"G1T12R_kegg_enrichment_results.txt",sep = "\t",quote = F)


# G2T12L kegg

DElist_low <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype2_leaf_treat1_VS_treat2/G2_L_low_expression_id.txt")
DElist_up <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype2_leaf_treat1_VS_treat2/G2_L_up_expression_id.txt")
colnames(DElist_up) <- "gene_id"
id <- rbind(DElist_low,DElist_up)
id$gene_id


x <- enricher(id$gene_id,
              TERM2GENE = pathway_gene,
              TERM2NAME = pathway_name,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.1)
barplot(x)
dotplot(x)
out <- as.data.frame(x)

write.table(out,"G2T12L_kegg_enrichment_results.txt",sep = "\t",quote = F)


# G2T12R kegg
DElist_low <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype2_root_treat1_VS_treat2/G2_R_low_expression_id.txt")
DElist_up <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype2_root_treat1_VS_treat2/G2_R_up_expression_id.txt")
colnames(DElist_up) <- "gene_id"
id <- rbind(DElist_low,DElist_up)
id$gene_id


x <- enricher(id$gene_id,
              TERM2GENE = pathway_gene,
              TERM2NAME = pathway_name,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.1)
barplot(x)
dotplot(x)
out <- as.data.frame(x)

write.table(out,"G2T12R_kegg_enrichment_results.txt",sep = "\t",quote = F)


# G1T12L GO

term2gene <- read.table("../ref/term_gene1.txt",header = F,sep = "\t")

term_name <- read.csv("../ref/term_name.txt",sep = ";",header = F)

term2class <- term_name[,c(1,3)]
term2name <- term_name[,1:2]
DElist_low <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype1_leaf_treat1_VS_treat2/G1_L_low_expression_id.txt")
DElist_up <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype1_leaf_treat1_VS_treat2/G1_L_up_expression_id.txt")
colnames(DElist_up) <- "gene_id"
id <- rbind(DElist_low,DElist_up)
id$gene_id


x <- enricher(id$gene_id,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.1)


out <- as.tibble(x)


out1 <- out %>% left_join(term2class,by = c("ID" = "V1"))

# 先按V3排序，再按Count排序，并创建一列用于排序
out1 <- out1 %>% arrange(V3,Count) %>% mutate(id = seq(1,nrow(out1)))


out1$Description <- factor(out1$Description,levels = out1$Description[order(out1$id)])

ggplot(data = out1,aes(x = Description,y = Count)) +
  geom_bar(stat = "identity",aes(fill = V3)) + 
  coord_flip() +
  theme_test() +
  guides(fill = guide_legend(title = "type",reverse = T)) +
  xlab("GO term")
 

# G1T12R go

DElist_low <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype1_root_treat1_VS_treat2/G1_R_low_expression_id.txt")
DElist_up <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype1_root_treat1_VS_treat2/G1_R_up_expression_id.txt")
colnames(DElist_up) <- "gene_id"
id <- rbind(DElist_low,DElist_up)
id$gene_id


x <- enricher(id$gene_id,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.1)


out <- as.tibble(x)


out1 <- out %>% left_join(term2class,by = c("ID" = "V1"))

# 先按V3排序，再按Count排序，并创建一列用于排序
out1 <- out1 %>% arrange(V3,Count) %>% mutate(id = seq(1,nrow(out1)))


out1$Description <- factor(out1$Description,levels = out1$Description[order(out1$id)])

ggplot(data = out1,aes(x = Description,y = Count)) +
  geom_bar(stat = "identity",aes(fill = V3)) + 
  coord_flip() +
  theme_test() +
  guides(fill = guide_legend(title = "type",reverse = T)) +
  xlab("GO term")

# G2T12L go

DElist_low <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype2_leaf_treat1_VS_treat2/G2_L_low_expression_id.txt")
DElist_up <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype2_leaf_treat1_VS_treat2/G2_L_up_expression_id.txt")
colnames(DElist_up) <- "gene_id"
id <- rbind(DElist_low,DElist_up)
id$gene_id


x <- enricher(id$gene_id,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.1)


out <- as.tibble(x)


out1 <- out %>% left_join(term2class,by = c("ID" = "V1"))

# 先按V3排序，再按Count排序，并创建一列用于排序
out1 <- out1 %>% arrange(V3,Count) %>% mutate(id = seq(1,nrow(out1)))


out1$Description <- factor(out1$Description,levels = out1$Description[order(out1$id)])

ggplot(data = out1,aes(x = Description,y = Count)) +
  geom_bar(stat = "identity",aes(fill = V3)) + 
  coord_flip() +
  theme_test() +
  guides(fill = guide_legend(title = "type",reverse = T)) +
  xlab("GO term")

# G2T12R go

DElist_low <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype2_root_treat1_VS_treat2/G2_R_low_expression_id.txt")
DElist_up <- read.table("../Results_to_fareeha/Different_analysis_results/Genetype2_root_treat1_VS_treat2/G2_R_up_expression_id.txt")
colnames(DElist_up) <- "gene_id"
id <- rbind(DElist_low,DElist_up)
id$gene_id


x <- enricher(id$gene_id,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.1)


out <- as.tibble(x)


out1 <- out %>% left_join(term2class,by = c("ID" = "V1"))

# 先按V3排序，再按Count排序，并创建一列用于排序
out1 <- out1 %>% arrange(V3,Count) %>% mutate(id = seq(1,nrow(out1)))


out1$Description <- factor(out1$Description,levels = out1$Description[order(out1$id)])

ggplot(data = out1,aes(x = Description,y = Count)) +
  geom_bar(stat = "identity",aes(fill = V3)) + 
  coord_flip() +
  theme_test() +
  guides(fill = guide_legend(title = "type",reverse = T)) +
  xlab("GO term")
