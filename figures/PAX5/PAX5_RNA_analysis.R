####################################################################################################
## Script purpose: RNA analysis on combined PAX5 repeat3
## Author: Zepeng Mu
## Date: Mon Jan  6 20:24:46 2025
####################################################################################################
library(tidyverse)
library(data.table)
library(Seurat)
library(SeuratObject)
library(matrixStats)
library(sparseMatrixStats)
library(SCpubr)
library(harmony)
"%&%" <- function(a, b) paste0(a, b)

basicProcess <- function(object = NULL, nfeatures = 1000) {
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeatures)
  
  allGenes <- rownames(object)
  object <- ScaleData(object, features = allGenes)
  
  object <- RunPCA(object, features = VariableFeatures(object = object))
  
  DimPlot(object, reduction = "pca") + NoLegend()
  
  object <- FindNeighbors(object, dims = 1:10)
  
  return(object)
}

# Load data ----
smps <- "PAX5_"%&%c("A1", "A2", "B1", "B2", "M1", "M2")

pax5Cmb <- read_rds("../data/pax5Cmb.rds")

geno <- fread("../data/allCells_haplo_geno_AB.txt.gz")

metaDf <- pax5Cmb@meta.data
metaDf <- metaDf %>% 
  dplyr::left_join(geno %>% dplyr::select(smp, well, A.T_84:A.T_108, B.T_64:cellType, -bc),
                   by = c("orig.ident" = "smp", "well"))

rownames(metaDf) <- colnames(pax5Cmb)
pax5Cmb@meta.data <- metaDf

# Keep cells that have genotype data
pax5Cmb <- pax5Cmb[, pax5Cmb$orig.ident%&%pax5Cmb$well %in% (geno$smp%&%geno$well)]

# RNA ----
Idents(pax5Cmb) <- "orig.ident"
VlnPlot(pax5Cmb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nFeature_ADT", "nCount_ADT"), ncol = 3)
pax5Cmb <- subset(pax5Cmb, subset = nCount_RNA < 15e4 & nCount_ADT < 25000)
VlnPlot(pax5Cmb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nFeature_ADT", "nCount_ADT"), ncol = 3)

pax5Cmb <- basicProcess(pax5Cmb, nfeatures = 1500)
pax5Cmb <- FindClusters(pax5Cmb, resolution = 0.2)
pax5Cmb <- RunUMAP(pax5Cmb, dims = 1:15, n.neighbors = 15, spread = 0.8, min.dist = 1)

## Harmony ----
pax5Cmb[["RNA"]] <- split(pax5Cmb[["RNA"]], f = pax5Cmb$orig.ident)
pax5Cmb <- IntegrateLayers(pax5Cmb, method = HarmonyIntegration, orig.reduction = "pca",
                           new.reduction = "harmony")

pax5Cmb <- FindNeighbors(pax5Cmb, reduction = "harmony", dims = 1:15)
pax5Cmb <- FindClusters(pax5Cmb, resolution = 0.2, cluster.name = "harmony_clusters")
pax5Cmb <- RunUMAP(pax5Cmb, reduction = "harmony", dims = 1:15, n.neighbors = 15,
                   reduction.name = "har_umap", spread = 0.8, min.dist = 0.6)

g1 <- do_DimPlot(pax5Cmb, pt.size = 0.1, plot.title = "Condition", reduction = "har_umap",
                 group.by = "cellType") +
  theme_zm(base_size = 10) +
  labs(x = "Harmony UMAP 1", y = "Harmony UMAP 2") +
  theme(aspect.ratio = 1,
        legend.position = "bottom")

g2 <- do_DimPlot(pax5Cmb, pt.size = 0.1, plot.title = "Harmony clusters", reduction = "har_umap",
                 group.by = "harmony_clusters") +
  theme_zm() +
  labs(x = "Harmony UMAP 1", y = "Harmony UMAP 2") +
  theme(aspect.ratio = 1,
        legend.position = "bottom")

g1 | g2

pax5Cmb <- JoinLayers(pax5Cmb)

write_rds(pax5Cmb, "../data/pax5Cmbmony.rds")

# Cluster analysis ----
DefaultAssay(pax5Cmb) <- "RNA"

pax5Cmb$batch <- str_sub(pax5Cmb$orig.ident, 6, 6)
pax5Cmb$condition <- pax5Cmb$batch %&% "-" %&% pax5Cmb$cellType

do_DimPlot(pax5Cmb,
           pt.size = 0.1, plot.title = "Condition", reduction = "har_umap",
           group.by = "cellType", colors.use = c("NTC" = "grey", "Edited" = "orange2")) +
  labs(x = "Harmony UMAP 1", y = "Harmony UMAP 2") +
  theme(aspect.ratio = 1, legend.position = "bottom")

do_DimPlot(pax5Cmb,
           pt.size = 0.1, plot.title = "Cluster", reduction = "har_umap",
           group.by = "harmony_clusters",
           colors.use = c("0" = "darkgreen", "1" = "darkred", "2" = "blue2")) +
  labs(x = "Harmony UMAP 1", y = "Harmony UMAP 2") +
  theme(aspect.ratio = 1, legend.position = "bottom")

pax5Cmb@meta.data %>%
  filter(cellType != "NTC") %>% 
  dplyr::group_by(harmony_clusters, batch) %>%
  tally() %>%
  ggplot(aes(x = harmony_clusters, y = n, fill = batch)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(
    limits = c("A", "B", "M"), values = c("yellow2", "green4", "blue2"),
    labels = c("A", "B", "A+B"), name = ""
  ) +
  labs(x = "Clusters", y = "Proportion") +
  theme(legend.position = "top")

## Markers ----
Idents(pax5Cmb) <- pax5Cmb$harmony_clusters

diffGene <- FindMarkers(pax5Cmb,
                        assay = "RNA", logfc.threshold = 0,
                        min.pct = 0.25, ident.1 = "0", ident.2 = "1")

diffGeneSig <- diffGene %>%
  filter(p_val_adj < 0.01 & abs(avg_log2FC) > 1)

diffGene1 <- FindMarkers(pax5Cmb,
                         assay = "RNA", logfc.threshold = 0,
                         min.pct = 0.25, ident.1 = "0", ident.2 = "2")

diffGeneSig1 <- diffGene1 %>%
  filter(p_val_adj < 0.01 & abs(avg_log2FC) > 1)

## CD79B FCER2 IGHG3 VPREB3
tmpGene <- "VPREB3" # Change accordingly
do_ViolinPlot(pax5Cmb,
              features = tmpGene, font.size = 10, pt.size = 0.1, line_width = 0, order = F,
              plot_boxplot = F, colors.use = c("0" = "darkgreen", "1" = "darkred", "2" = "blue2")) +
  labs(x = "Cluster", y = str_glue("*{tmpGene}*")) +
  theme(legend.position = "bottom", axis.title.y = ggtext::element_markdown())
