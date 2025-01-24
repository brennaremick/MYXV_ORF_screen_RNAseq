library(plyr)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)

# read in data
sample_info <- read.table("barcodes_samples.txt", header=TRUE) # has barcode info and sample names
combined_df <- read.table("screen_raw_reads.txt", header=TRUE, row.names = 1) # raw counts for all samples

# iterate through each ORF and make a PCA plot comparing ORF and mCherry

# identify all conditions
levels <- unique(sample_info$condition)

# remove samples that do not have at least 2 biological replicates ("single" = samples that only have 1 replicate)
# also remove "mCherry" from levels
remove <- c("single", "mCherry")
levels <- levels[!unlist(levels) %in% unlist(remove)]

# remove MYXV genes from combined_df
combined_df_human <- combined_df[grepl("ENSG", rownames(combined_df)), ]

# iterate through all ORFs (levels) and create a PCA plot compared to mCherry
pdf("pca_plots_myxv_screen.pdf")
for(i in 1:length(levels)) {
  info <- sample_info[ which (sample_info$condition == "mCherry" | sample_info$condition == levels[i]), ]
  subset <- subset(combined_df_human, select = info$sample)
  dds <- DESeqDataSetFromMatrix(countData = subset, colData = info, design = ~condition)
  keep_rows <- rowSums(counts(dds)) > 1
  dds <- dds[keep_rows,]
  vsd <- vst(dds, blind = FALSE)
  print(plotPCA(vsd, intgroup = c("condition")))
}
dev.off()


# iterate through all the ORFs and perform differential gene expression analysis compared to mCherry
df_list <- list()
for(i in 1:length(levels)) {
  info <- sample_info[ which (sample_info$condition == "mCherry" | sample_info$condition == levels[i]), ]
  subset <- subset(combined_df_human, select = info$sample)
  dds <- DESeqDataSetFromMatrix(countData = subset, colData = info, design = ~condition)
  keep_rows <- rowSums(counts(dds)) > 0
  dds <- dds[keep_rows,]
  dds$condition <- relevel(dds$condition, ref = "mCherry")
  dds <-DESeq(dds)
  res <- as.data.frame(results(dds))
  res$comparison <-levels[i] # add a column with the ORF name
  res_ordered <- res[order(-res$log2FoldChange), ]
  df_list[[i]] <- res_ordered   # Append the dataframe to the list
}

df_dge <- do.call(rbind, lapply(df_list, as.data.frame)) # convert df_list into a dataframe

# add symbol column to dataframe with gene symbols
ens.str <- substr(rownames(df_dge), 1, 15)
df_dge$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")


