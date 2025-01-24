library(DESeq2)
library(EnhancedVolcano)
library(org.Hs.eg.db)

# identify differentially expressed genes in M003.1-expressing samples compared to mCherry

# read in data
sample_info <- read.table("barcodes_samples.txt", header=TRUE) # has barcode info and sample names
combined_df <- read.table("screen_raw_reads.txt", header=TRUE, row.names = 1) # raw counts for all samples

# remove MYXV genes from combined_df
combined_df_human <- combined_df[grepl("ENSG", rownames(combined_df)), ]

# subset just M003.1 and mCherry samples
info <- sample_info[ which (sample_info$condition == "mCherry" | sample_info$condition == "m003_1L"), ]
subset <- subset(combined_df_human, select = info$sample)

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = subset, colData = info, design = ~condition)

# pre-filtering to remove low count genes
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]

# choose reference level; compare M003.1 samples to mCherry
dds$condition <- relevel(dds$condition, ref = "mCherry")

dds <- DESeq(dds)
res <- results(dds)
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

# make volcano plot (Fig S3)
EnhancedVolcano(res,
                lab = res$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 2,
                pCutoff = 10e-5,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                gridlines.minor = FALSE,
                gridlines.major = FALSE,
                max.overlaps = 9)