library(plyr)
library(dplyr)

# use this script to take the zUMIs output and create a dataframe that is compatible with DESeq2

# import all rds zUMIs output files
a = readRDS("BR_myxv_human_A.dgecounts.rds")
b = readRDS("BR_myxv_human_B.dgecounts.rds") 
c = readRDS("BR_myxv_human_C.dgecounts.rds") 
d = readRDS("BR_myxv_human_D.dgecounts.rds") 
e = readRDS("BR_myxv_human_E.dgecounts.rds") 
f = readRDS("BR_myxv_human_F.dgecounts.rds") 
g = readRDS("BR_myxv_human_G.dgecounts.rds") 
h= readRDS("BR_myxv_human_H.dgecounts.rds") 
j= readRDS("BR_myxv_human_J.dgecounts.rds") 
k= readRDS("BR_myxv_human_K.dgecounts.rds") 
l= readRDS("BR_myxv_human_L.dgecounts.rds") 
m= readRDS("BR_myxv_human_M.dgecounts.rds") 
n= readRDS("BR_myxv_human_N.dgecounts.rds") 
o= readRDS("BR_myxv_human_O.dgecounts.rds") 
p= readRDS("BR_myxv_human_P.dgecounts.rds") 
q= readRDS("BR_myxv_human_Q.dgecounts.rds") 
r= readRDS("BR_myxv_human_R.dgecounts.rds") 
s= readRDS("BR_myxv_human_S.dgecounts.rds") 
t= readRDS("BR_myxv_human_T.dgecounts.rds") 

rds_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t")

for(i in 1:length(rds_names)) {
  assign(paste0("umi_inex_all_", rds_names[i]),as.data.frame(as.matrix(eval(as.symbol(rds_names[i]))$umicount$inex$all)))
}

# tagmentation was performed with 4-5 technical replicates per library
# add reads for each sample across tagmentation technical replicates (use this function for 4 tech rep)
add_tech_replicates_4 <- function(df1, df2, df3, df4) {
  pp <- cbind(names=c(rownames(df1), rownames(df2), rownames(df3), rownames(df4)), 
              rbind.fill(list(df1, df2, df3, df4)))
  ddply(pp, .(names), function(x) colSums(x[,-1], na.rm = TRUE))
}

# add reads for each sample across technical replicates (use this function for 5 tech rep)
add_tech_replicates_5 <- function(df1, df2, df3, df4, df5) {
  pp <- cbind(names=c(rownames(df1), rownames(df2), rownames(df3), rownames(df4),rownames(df5)), 
              rbind.fill(list(df1, df2, df3, df4, df5)))
  ddply(pp, .(names), function(x) colSums(x[,-1], na.rm = TRUE))
}

dup1_p1 <- add_tech_replicates_4(umi_inex_all_a, umi_inex_all_e, umi_inex_all_m, umi_inex_all_q) # dup1_p1 = plate 1 (96-well plate); biological replicate 1
dup1_p2 <- add_tech_replicates_5(umi_inex_all_b, umi_inex_all_f, umi_inex_all_j, umi_inex_all_n, umi_inex_all_r) # dup1_p2 = plate 2 (96-well plate); biological replicate 1
dup2_p1 <- add_tech_replicates_5(umi_inex_all_c, umi_inex_all_g, umi_inex_all_k, umi_inex_all_o, umi_inex_all_s) # dup2_p1 = plate 1 (96-well plate); biological replicate 2
dup2_p2 <- add_tech_replicates_5(umi_inex_all_d, umi_inex_all_h, umi_inex_all_l, umi_inex_all_p, umi_inex_all_t) # dup2_p2 = plate 2 (96-well plate); biological replicate 2

# add a unique label to each sample name
colnames(dup1_p1) <- paste(colnames(dup1_p1),"1",sep="_")
colnames(dup1_p2) <- paste(colnames(dup1_p2),"2",sep="_")
colnames(dup2_p1) <- paste(colnames(dup2_p1),"3",sep="_")
colnames(dup2_p2) <- paste(colnames(dup2_p2),"4",sep="_")

# assign row names and remove first column
rownames(dup1_p1) <- dup1_p1[,1] #Assigning row names from 1st column 
dup1_p1[,1] <- NULL #Removing the first column
rownames(dup1_p2) <- dup1_p2[,1]
dup1_p2[,1] <- NULL 
rownames(dup2_p1) <- dup2_p1[,1] 
dup2_p1[,1] <- NULL 
rownames(dup2_p2) <- dup2_p2[,1] 
dup2_p2[,1] <- NULL

# Merge into 1 dataframe
combined_df <- merge(dup1_p1, dup1_p2, by="row.names",all=TRUE)
rownames(combined_df) <- combined_df[,1]  
combined_df[,1] <- NULL
combined_df <- merge(combined_df, dup2_p1, by="row.names",all=TRUE)
rownames(combined_df) <- combined_df[,1] 
combined_df[,1] <- NULL
combined_df <- merge(combined_df, dup2_p2, by="row.names",all=TRUE)
rownames(combined_df) <- combined_df[,1] 
combined_df[,1] <- NULL

# convert all NA to 0
combined_df[is.na(combined_df)] <- 0

