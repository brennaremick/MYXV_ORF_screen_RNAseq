rm(list = ls())
library(tidyverse)
library(dplyr)

# read in the data
# DESeq2 was used to perform differential gene expression analysis for each ORF compared to mCherry (see analysis_screen.R script)
dat = read.table("MYXV_screen_DESeq2.txt", header=TRUE, row.names = 1)

# remove conditions where the cells were dying/very unhealthy looking; it is likely the transduction didn't work
sick <- c("m002L",
          "m006L",
          "m008_1L",
          "m008L",
          "m009L",
          "m028L",
          "m032R",
          "m034L",
          "m044R",
          "m045L",
          "m046L",
          "m054R",
          "m047R",
          "m068R",
          "m072L",
          "m088L",
          "m099L",
          "m111R",
          "m127L",
          "m133R",
          "m144R")

dat <- dat[!(dat$comparison %in% sick), ]

# replace NA values with 1
ggdf = dat %>% replace_na(list(pvalue = 1, padj = 1)) #%>%
# mutate(padj = max(padj, 1e-200))

ggdf$padj[which(!is.finite(-log10(ggdf$padj)))]

# Select the top 100 rows with the smallest 'padj' value within each group
ggdf_subset = ggdf %>% group_by(comparison) %>% slice_min(padj, n = 100) %>% ungroup()

# Plot the empirical cumulative distribution (ECDF) of the absolute log2FC for the top 100 genes with the lowest 'padj' values for each condition 
ggplot(ggdf_subset, aes(x = abs(log2FoldChange), color = comparison)) +
  stat_ecdf() + scale_color_discrete(guide = NULL)



