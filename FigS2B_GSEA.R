library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)

# read in data
sample_info <- read.table("barcodes_samples.txt", header=TRUE) # has barcode info and sample names
combined_df <- read.table("screen_raw_reads.txt", header=TRUE, row.names = 1) # raw counts for all samples

# remove MYXV genes from combined_df
combined_df_human <- combined_df[grepl("ENSG", rownames(combined_df)), ]

# subset each interesting ORF and perform GSEA
# example:
info <- sample_info[ which (sample_info$condition == "mCherry" | sample_info$condition == "m003_1L"), ]
subset <- subset(combined_df_human, select = info$sample)
dds <- DESeqDataSetFromMatrix(countData = subset, colData = info, design = ~condition)
keep <- rowSums(counts(dds)) >= 10 # pre-filtering to remove low count genes
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "mCherry") # choose reference level; comapre to  mCherry
dds <- DESeq(dds)
res <- results(dds)
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

mydata.df.sub <- data.frame(geneID = res$symbol, LogFC = res$log2FoldChange)

# construct a named vector
mydata.gsea <- mydata.df.sub$LogFC 
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res.H <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_H, verbose=FALSE)
myGSEA.df.H <- as_tibble(myGSEA.res.H@result)
myGSEA.df.H.top5 <- myGSEA.df.H %>% 
  top_n(5, NES) # identify top 5 enriched pathways

# plot top 5 enriched pathways; use the sample scale for all plots
plot_M003_1 <- ggplot(myGSEA.df.H.top5, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(
    colors = c("red", "blue"), 
    trans = "log10",                    
    limits = c(5e-010, 0.05),            # Define scale limits
    breaks = c(5e-08, 5e-06, 5e-04), # Define scale breakpoints
    labels = c("5e-08", "5e-06","5e-04") # Define scale labels
  ) +
  ylim(limits = c(-3, 3)) +
  labs(y = "Normalized enrichment score", x = "", fill = "p.adjust") +
  coord_flip() +
  guides(fill = guide_colorbar(reverse = TRUE)) +
  theme_bw() 