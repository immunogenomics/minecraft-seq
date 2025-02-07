####################################################################################################
## Script purpose: Bulk RPL8 analysis
## Author: Zepeng Mu
## Date: Mon Jan 27 15:06:10 2025
####################################################################################################
library(data.table)
library(tidyverse)
library(matrixStats)
library(Matrix)
library(DESeq2)
library(ggtext)
library(ggrepel)
library(qvalue)
library(tidyplots)
library(ggrastr)
"%&%" <- function(a, b) paste0(a, b)

## RNA-seq count data ----
cnt <- fread("../data/gene_name_counts.txt")
colnames(cnt) <- str_remove_all(colnames(cnt), "/data/srlab1/zmu/revision/data/RPL8/bulk_RNA/|/Aligned.sortedByCoord.out.bam")
cnt <- cnt %>% filter(!startsWith(Geneid, fixed("MT-")))

cntMtrx <- as.matrix(cnt %>% dplyr::select(-(Geneid:Length)))
rownames(cntMtrx) <- cnt$Geneid

normCntMtrx <- t(t(cntMtrx) / colSums2(cntMtrx) * 1e6)
sum(rowMeans2(normCntMtrx) > 0.5)

gene2keep <- names(which(rowSums2(cntMtrx[, "RPL8_E"%&%1:5] > 0) >= 3 &
                           rowSums2(cntMtrx[, "RPL8_N"%&%1:5] > 0) >= 3))

cntFilter <- cntMtrx[intersect(gene2keep, names(which(rowMeans2(normCntMtrx) >= 1))), ]
dim(cntFilter)

metaDt <- data.frame(name = colnames(cntFilter)) %>% 
  separate(name, c("exp", "condSmp"), sep = "_", remove = F) %>% 
  mutate(cond = factor(str_remove(condSmp, "[0-9]"), levels = c("N", "E")),
         smp = as.factor(str_sub(condSmp, -1)))

# Diff analysis ----
dds <- DESeqDataSetFromMatrix(countData = cntFilter, colData = metaDt, design = ~ cond)
dds <- DESeq(dds)

res <- results(dds, name = "cond_E_vs_N")
resDf <- as.data.frame(res)
resDf$gene <- rownames(resDf)
resDf$qval <- qvalue(resDf$pvalue)$qvalues

res1 <- lfcShrink(dds, coef = "cond_E_vs_N", type = "ashr")
resDf1 <- as.data.frame(res1)
resDf1$gene <- rownames(resDf1)
resDf1$qval <- qvalue(resDf1$pvalue)$qvalues

fwrite(resDf1, sep = "\t", row.names = F, col.names = T, quote = F,
       file = "../results/diff_bulk_deseq2_ashr.txt.gz")

gene2plot <- c("RGS1", "PFKFB4", "FGF12", "EGR1", "EGR2", "MS4A7", "MS4A14", "CD69", "DHRS2")

resDf1 %>% 
  mutate(sig = case_when(log2FoldChange > 1 & qval < 0.05 ~ "red2",
                         log2FoldChange < -1 & qval < 0.05 ~ "blue2",
                         T ~ "grey")) %>% 
  ggplot(aes(log2FoldChange, -log10(pvalue))) +
  theme_zm() +
  geom_point_rast(aes(col = sig), size = 0.3, raster.dpi = 960, scale = 0.8) +
  geom_text_repel(data = res %>% filter(gene %in% gene2plot), aes(label = gene),
                  fontface = "italic", size = 2.5) +
  geom_hline(yintercept = -log10(pvalCutoff), lty = 2) +
  geom_vline(xintercept = c(-1, 1), lty = 2) +
  scale_color_identity() +
  labs(x = "log<sub>2</sub>FC", y = "-log<sub>10</sub>(P)") +
  theme(aspect.ratio = 1,
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown())

table(resDf$log2FoldChange > 0.5 & resDf$qval < 0.05)
table(resDf$log2FoldChange < -0.5 & resDf$qval < 0.05)

metaDt %>% 
  mutate(tmpCnt = log1p(normCntMtrx["RPL8", .data$name])) %>% 
  tidyplot(x = cond, y = tmpCnt) %>% 
  adjust_font(fontsize = 12) %>% 
  add_boxplot(box_width = 0.4,show_outliers = F) %>% 
  add_data_points(size = 2) %>% 
  adjust_x_axis_title("") %>% 
  adjust_y_axis("RPL8 log(CPM+1)") %>% 
  adjust_x_axis(labels = c("NTC", "Edited"))
