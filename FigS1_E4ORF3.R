library(DESeq2)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)

# identify differential expressed genes in E4ORF3-expressing samples compared to mCherry

# read in data
sample_info <- read.table("barcodes_samples.txt", header=TRUE) # has barcode info and sample names
combined_df <- read.table("screen_raw_reads.txt", header=TRUE, row.names = 1) # raw counts for all samples

# remove MYXV genes from combined_df
combined_df_human <- combined_df[grepl("ENSG", rownames(combined_df)), ]

# subset just E4ORF and mCherry samples
info <- sample_info[ which (sample_info$condition == "mCherry" | sample_info$condition == "E4ORF3"), ]
subset <- subset(combined_df_human, select = info$sample)

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = subset, colData = info, design = ~condition)

# pre-filtering to remove low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# choose reference level; compare E4ORF3 samples to mCherry
dds$condition <- relevel(dds$condition, ref = "mCherry")

dds <- DESeq(dds)
res <- results(dds)
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

# make volcano plot (Fig S1A)
EnhancedVolcano(res,
                lab = res$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 2,
                pCutoff = 10e-10,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                gridlines.minor = FALSE,
                gridlines.major = FALSE,
                max.overlaps = 9)

# GSEA (Fig S1B)
# Download Hallmark gene sets from MSigDB
hs_gsea_H <- msigdbr(species = "Homo sapiens", 
                     category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

# Pull out just the columns corresponding to gene symbols and LogFC 
mydata.df.sub <- data.frame(geneID = res$symbol, LogFC = res$log2FoldChange)

# Construct a named vector
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res.H <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_H, verbose=FALSE)
myGSEA.df.H <- as_tibble(myGSEA.res.H@result)

# plot
ggplot(myGSEA.df.H, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(y = "Normalized enrichment score", x="") +
  coord_flip() + # Flip the coordinates to make it horizontal
  guides(fill = guide_colorbar(reverse = TRUE)) +
  theme_bw()

