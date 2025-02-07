####################################################################################################
## Script purpose: QC and data analysis on bulk RNA-seq data
## Author: Zepeng Mu
## Date: Mon Nov 11 11:18:41 2024
####################################################################################################
library(matrixStats)
library(data.table)
library(tidyverse)
library(Matrix)
library(DESeq2)
library(ggtext)
library(ggrepel)
library(qvalue)
library(fgsea)
library(tidyplots)
library(ggrastr)
"%&%" <- function(a, b) paste0(a, b)

# RNA-seq ----
## Load cell-by-gene count matrix
cnt <- fread("/data/srlab1/zmu/revision/data/CD45KO_R/bulk/RNA/featureCounts.txt")
colnames(cnt) <- str_remove_all(colnames(cnt), "/data/srlab1/zmu/revision/data/CD45KO_R/bulk/RNA/|/STARsolo/Aligned.sortedByCoord.out.bam")
colnames(cnt)[colnames(cnt) == "CD45_KO4"] <- "CD45_NTC4" # Fix label switch
cnt <- cnt %>% filter(!startsWith(Geneid, fixed("MT-"))) # Remove mitochondria gene

# Remove low quality library
cntMtrx <- as.matrix(cnt %>% dplyr::select(-(Geneid:Length), -CD45_KO6))
rownames(cntMtrx) <- cnt$Geneid

normCntMtrx <- t(t(cntMtrx) / colSums2(cntMtrx) * 1e6)
sum(rowMeans2(normCntMtrx) > 0.5)

# Remove lowly expressed genes
gene2keep <- names(which(rowSums2(cntMtrx[, c("CD45_KO1", "CD45_KO2", "CD45_KO3", "CD45_KO5")] > 0) >= 3 &
                           rowSums2(cntMtrx[, "CD45_NTC"%&%2:6] > 0) >= 3))

cntFilter <- cntMtrx[intersect(gene2keep, names(which(rowMeans2(normCntMtrx) >= 1))), ]
dim(cntFilter)

# Create metadata
metaDt <- data.frame(name = colnames(cntFilter)) %>% 
  separate(name, c("exp", "condDonor"), sep = "_", remove = F) %>% 
  mutate(cond = factor(str_remove(condDonor, "[0-9]"), levels = c("NTC", "KO")),
         donor = as.factor(str_sub(condDonor, -1)))

# Diff analysis ----
dds <- DESeqDataSetFromMatrix(countData = cntFilter, colData = metaDt, design = ~ cond + donor)
dds <- DESeq(dds, test = "LRT", reduced = ~ donor)

res <- results(dds, name = "cond_KO_vs_NTC")
resDf <- as.data.frame(res)
resDf$gene <- rownames(resDf)
resDf$qval <- qvalue(resDf$pvalue)$qvalues

fwrite(resDf, sep = "\t", row.names = F, col.names = T, quote = F,
       file = "../results/diff_bulk_deseq2.txt.gz")

sum(resDf$log2FoldChange > 0.5 & resDf$qval < 0.05)
sum(resDf$log2FoldChange < -0.5 & resDf$qval < 0.05)

sum(resDf$log2FoldChange > 1 & resDf$qval < 0.05)
sum(resDf$log2FoldChange < -1 & resDf$qval < 0.05)

gene2plot <- c("PTPRC", "CD4", "TNFRSF4", "TRBC2", "ANXA3", "CHCHD10", "DCTN3",
               "IFI6", "ISG20", "RPS3", "ISG15", "RPL13", "PSMB8")

resDf %>% 
  mutate(sig = case_when(log2FoldChange > 0.5 & qval < 0.05 ~ "red2",
                         log2FoldChange < -0.5 & qval < 0.05 ~ "blue2",
                         T ~ "grey")) %>% 
  ggplot(aes(log2FoldChange, -log10(pvalue))) +
  geom_point_rast(aes(col = sig), size = 0.2, raster.dpi = 960, scale = 0.5) +
  geom_text_repel(data = resDf %>% filter(gene %in% gene2plot & log2FoldChange > 0),
                  aes(label = gene), size = 2, xlim = c(2, Inf), force = 3, fontface = "italic", segment.size = 0.2) +
  geom_text_repel(data = resDf %>% filter(gene %in% gene2plot & log2FoldChange < 0),
                  aes(label = gene), size = 2, xlim = c(-Inf, -2), force = 3, fontface = "italic", segment.size = 0.2) +
  geom_hline(yintercept = -log10(max(resDf$pvalue[resDf$qval < 0.05])), lty = 2) +
  geom_vline(xintercept = c(-0.5, 0.5), lty = 2) +
  scale_color_identity() +
  labs(x = "log<sub>2</sub>FC", y = "-log<sub>10</sub>(P)") +
  theme(aspect.ratio = 1,
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown())

metaDt %>% 
  mutate(tmpCnt = log1p(normCntMtrx["PTPRC", .data$name])) %>% 
  tidyplot(x = cond, y = tmpCnt) %>% 
  adjust_font(fontsize = 12) %>% 
  add_boxplot(box_width = 0.4,show_outliers = F) %>% 
  add_data_points(size = 2) %>% 
  adjust_x_axis_title("") %>% 
  adjust_y_axis("PTPRC log(CPM+1)")

## GSEA analysis ----
c7gmt <- gmtPathways("../data/immsigDB_C7.gmt")
sigRank <- -log10(resDf$pvalue) * sign(resDf$log2FoldChange)
names(sigRank) <- resDf$gene
sigRank <- sort(sigRank)
fgseaRes <- fgsea(c7gmt, sigRank, minSize = 15, maxSize = 500)

# Format fGSEA output
fgseaRes1 <- fgseaRes %>% 
  filter(padj < 0.1) %>% 
  rowwise() %>% 
  mutate(leadingEdge = str_flatten(unlist(leadingEdge), collapse = ",")) %>% 
  ungroup()

fgseaRes1 %>% 
  dplyr::group_by(sign(NES)) %>% 
  slice_max(abs(NES), n = 10) %>% 
  ungroup() %>% 
  arrange(NES) %>% 
  mutate(pathway = str_to_title(str_replace_all(pathway, fixed("_"), " ")),
         pathway = factor(pathway, levels = pathway)) %>% 
  ggplot(aes(x = NES, y = pathway, fill = -log10(padj))) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  scale_fill_gradient(low = "darkblue", high = "red") +
  labs(y = "")

## GSEA plots
plotEnrichment(c7gmt[["GSE28726_NAIVE_VS_ACTIVATED_CD4_TCELL_UP"]], sigRank) +
  labs(title="NAIVE_VS_ACTIVATED_CD4_TCELL_UP (padj=3.02e-4)") +
  theme(axis.line = element_blank(), panel.grid.major = element_line(color = "grey90"))

plotEnrichment(c7gmt[["GSE3982_MEMORY_CD4_TCELL_VS_TH1_DN"]], sigRank) +
  labs(title="MEMORY_CD4_TCELL_VS_TH1_DN (padj=0.0054)") +
  theme(axis.line = element_blank(), panel.grid.major = element_line(color = "grey90"))
