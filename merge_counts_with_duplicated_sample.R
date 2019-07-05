

# 整理fareeha表达量 ----------------------------------------------------------


rm (list = ls())
setwd("/home/chenzhi/fareeha/count/results/")

file <- dir()
li <- vector("list",length = length(file))
names(li) <- file

count <- lapply(names(li),FUN = function(x) li[[x]] <- read.table(x,header = T,sep = "\t"))

dim(count[[1]])

library(tidyverse)

count1 <- Reduce(left_join,count)

name <- c("G1T1R","G1T2R",
          "G1T1L","G1T2L",
          "G2T1R","G2T2R",
          "G2T1L","G2T2L")
reli <- vector("list",length = length(name))
for (i in name) {
  reli[[i]] <- count1 %>% select(1,contains(i)) 
}

# 不同处理的基因count data.frame
G1T12R <- left_join(reli[["G1T1R"]],reli[["G1T2R"]])
G1T12L <- left_join(reli[["G1T1L"]],reli[["G1T2L"]])

G2T12R <- left_join(reli[["G2T1R"]],reli[["G2T2R"]])
G2T12L <- left_join(reli[["G2T1L"]],reli[["G2T2L"]])

#重命名
names(G1T12R)[2:7] <- c("G1T1R1","G1T1R2","G1T1R3",
                        "G1T2R1","G1T2R2","G1T2R3")

names(G1T12L)[2:7] <- c("G1T1L1","G1T1L2","G1T1L3",
                        "G1T2L1","G1T2L2","G1T2L3")

names(G2T12R)[2:7] <- c("G2T1R1","G2T1R2","G2T1R3",
                        "G2T2R1","G2T2R2","G2T2R3")

names(G2T12L)[2:7] <- c("G2T1L1","G2T1L2","G2T1L3",
                        "G2T2L1","G2T2L2","G2T2L3")


all_counts <- left_join(G1T12R,G1T12L)
all_counts <- left_join(all_counts,G2T12R)
all_counts <- left_join(all_counts,G2T12L)

write.table(all_counts,file = "all_counts.txt",sep = "\t",quote = F,
            row.names = F)

write.table(G1T12R,file = "G1T12R_counts.txt",sep = "\t",quote = F,
            row.names = F)

write.table(G1T12L,file = "G1T12L_counts.txt",sep = "\t",quote = F,
            row.names = F)
write.table(G2T12R,file = "G2T12R_counts.txt",sep = "\t",quote = F,
            row.names = F)
write.table(G2T12L,file = "G2T12L_counts.txt",sep = "\t",quote = F,
            row.names = F)
