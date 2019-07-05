
##汇总差异基因柱形图
setwd("../../../scripts/plot_script/")
rm(list = ls())
setwd("/home/chenzhi/fareeha/Results_to_fareeha/Different_analysis_results/Genetype1_leaf_treat1_VS_treat2/")

# G1 R+L 
dfL <- read.table("G1_L_low_expression_id.txt",header = T,row.names = 1)
 d1 <- nrow(dfL)

dfU <- read.table("G1_L_up_expression_id.txt",header = T,row.names = 1)
 d2 <- nrow(dfU)


dfL1 <- read.table("../Genetype1_root_treat1_VS_treat2/G1_R_low_expression_id.txt",header = T,row.names = 1)
 d3 <- nrow(dfL1)

dfU1 <- read.table("../Genetype1_root_treat1_VS_treat2/G1_R_up_expression_id.txt",header = T,row.names = 1)
 d4 <- nrow(dfU1)

# G2 R+L

dfL2 <- read.table("../Genetype2_leaf_treat1_VS_treat2/G2_L_low_expression_id.txt",header = T,row.names = 1)
dfL2
  d5 <- nrow(dfL2)
  
dfU2 <- read.table("../Genetype2_leaf_treat1_VS_treat2/G2_L_up_expression_id.txt",header = T,row.names = 1)
dfU2
  d6 <- nrow(dfU2)  

dfL3 <- read.table("../Genetype2_root_treat1_VS_treat2/G2_R_low_expression_id.txt",header = T,row.names = 1)
dfL3
  d7 <- nrow(dfL3)

dfU3 <- read.table("../Genetype2_root_treat1_VS_treat2/G2_R_up_expression_id.txt",header = T,row.names = 1)
dfU3
  d8 <- nrow(dfU3)

library(tidyverse)  

df <- tibble::tribble(~genetype,~class,~sample,~count,
                      "Genetype1","Up","Leaf",d2,
                      "Genetype1","Low","Leaf",d1,
                      "Genetype1","Up","Root",d4,
                      "Genetype1","Low","Root",d3,
                      "Genetype2","Up","Leaf",d6,
                      "Genetype2","Low","Leaf",d5,
                      "Genetype2","Up","Root",d8,
                      "Genetype2","Low","Root",d7,) 
df$class <- factor(df$class,levels = c("Up","Low"))
install.packages("ggsci")
library(ggsci)
library(RColorBrewer)
cols <- brewer.pal(1,"Set1")
cols1 <- brewer.pal(2,"Set1")
display.brewer.all()
df %>% ggplot(mapping = aes(x = sample,y = count,fill = class)) +
  geom_bar(stat = "identity",position = "dodge") + 
  facet_wrap(~genetype) + 
  theme_bw() + 
  scale_fill_aaas() + 
  scale_fill_discrete(breaks=c("Up","Low"))
  

