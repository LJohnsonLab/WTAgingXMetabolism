########################################################################
# Name: CCL_ST_MCC_SpaNorm_RPCA_WTAgingXMet_11.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 14. subsetting, subclustering, and annotating immune-vascular populations  
# Input Files: 
# Final Output Files: 
# Date created: 8/13/26
# Last updated: 8/21/26
# Author: Chloe Lucido
########################################################################

# Auto-install missing packages
required_packages <- c("Seurat", "future", "ggplot2", "readr", "dplyr", 
                       "data.table","tidyr", "tibble", "stringr", "forcats", "lubridate",  "readxl", "patchwork", 
                       "qs2", "RColorBrewer", "Polychrome", "purrr",
                       "presto", "glmGamPoi", "SCP", "pheatmap")

missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}


# Load Libraries ----
library(Seurat)
library(future)
library(ggplot2)
library(presto)
library(glmGamPoi)
library(readr)
library(dplyr)
library(data.table)
library(readxl)
library(patchwork)
library(qs2) 
library(RColorBrewer)
library(Polychrome)
library(purrr)
library(readxl)
library(tidyr)
library(tibble)
library(stringr)
library(speckle) # for propeller (cell prop stats)

# ONLY TO BE RUN ON MCC 

# Use the cores you requested
plan("multicore", workers = 8)
options(future.globals.maxSize = 400 * 1024^3)  # slightly under your --mem

# set seed
set.seed(42)

# to prevent this error: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so.
options(future.rng.onMisuse = "resolve")
options(future.seed = TRUE)


######### Paths ############
input_obj <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260617_SpaNorm_RPCA_scored.qs2" 
output_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/"
figs_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/SpaNorm_RPCA_figs/astro_only/"
cluster_markers <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/SpaNorm_cluster_markers/"
genes_to_remove_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/Genes_removed_from_astros.xlsx"
date <- "20260816_"
###########################



# 01. read in obj ----
SpaNorm_RPCA.obj <- qs_read(input_obj)

# 02. subset to astros ----
immune_vasc.obj <- subset(SpaNorm_RPCA.obj, subset = gen_celltype == c("Astrocyte", "Microglia", "T cell", "Vascular mural", "Pericyte", "Endothelial cell"))

# 03. save immune_vasc.obj to load locally ----
qs_save(immune_vasc.obj, paste0(output_path, date, "immune_vasc__SpaNorm_RPCA.qs2"))


# SANITY CHECK
dim(immune_vasc.obj)
#[1]   480 79833

table(immune_vasc.obj$gen_celltype)
#             Neuron        Astrocyte  Oligodendrocyte              OPC
#                  0            32281                0                0
#          Microglia           T cell Endothelial cell         Pericyte
#               9021              185            17738             4909
#     Vascular mural   Ependymal cell   Choroid plexus
#              15699                0                0

UMAP_1 <- DimPlot(immune_vasc.obj)

ggsave(filename = paste0(date, "immune_vasc_notprocessed_UMAP.png"), plot = UMAP_1, path = figs_path) 

# 04. process ----
# check default assay and adjust as needed 
DefaultAssay(immune_vasc.obj) <- "Xenium"

## 04a. normalize
immune_vasc.obj <- NormalizeData(immune_vasc.obj)

## 04b. find variable features 
immune_vasc.obj <- FindVariableFeatures(immune_vasc.obj)

## 04c. scale 
immune_vasc.obj <- ScaleData(immune_vasc.obj)

## 04d. PCA
immune_vasc.obj <- RunPCA(immune_vasc.obj, npcs = 50, features = rownames(immune_vasc.obj))

PCA_age <- PCAPlot(immune_vasc.obj, group.by = "Age")

ggsave(filename = paste0(date, "immune_vasc_PCA_byage.png"), plot = PCA_age, path = figs_path) 

# elbow plot to determine number of PCs to use
elbow_plot <- ElbowPlot(immune_vasc.obj, ndims = 50, reduction = "pca")

ggsave(filename = paste0(date, "immune_vasc_elbowplot.png"), plot = elbow_plot, path = figs_path) 

# save intermediate obj ----
qs_save(immune_vasc.obj, paste0(output_path, date, "SpaNorm_RPCA_immune_vasc_nodimreduc.qs2"))

# read in intermediate obj ----
immune_vasc.obj <- qs_read(paste0(output_path, "20260813_SpaNorm_RPCA_immune_vasc_nodimreduc.qs2"))

## 03e. RunUMAP ----
immune_vasc.obj <- RunUMAP(immune_vasc.obj, dims = 1:30)  

## 03f. FindNeighbors ----
immune_vasc.obj <- FindNeighbors(immune_vasc.obj, reduction = "pca", dims = 1:30)

## 03g. FindClusters ----
immune_vasc.obj <- FindClusters(immune_vasc.obj, resolution = 0.4, graph.name = "Xenium_snn")  

# save clusters to metadata ----
Idents(immune_vasc.obj) <- "seurat_clusters"

immune_vasc.obj$res04_clus <- Idents(immune_vasc.obj)

## 03h. Visualize UMAP with clusters ----
UMAP_res04 <- DimPlot(immune_vasc.obj, label = T, cols = "polychrome")

imgdim_splitclusters_res04 <- ImageDimPlot(immune_vasc.obj, split.by = "seurat_clusters", fov = "fov")

imgdim_res04 <- ImageDimPlot(immune_vasc.obj, group.by = "seurat_clusters", cols = "polychrome")

# join layers before finding markers
immune_vasc.obj <- JoinLayers(immune_vasc.obj)

## 03i. find markers ----
immunevasc.markers <- FindAllMarkers(immune_vasc.obj, group.by = "seurat_clusters", only.pos = T, logfc.threshold = 0.25)

immunevasc.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_immunevasc

write_csv(top5_immunevasc, paste0(cluster_markers, date, "SpaNorm_RPCA_immunevasc_subclusters_res04.csv"))  

# SANITY CHECK
table(immune_vasc.obj$res04_clus)
#     0     1     2     3     4     5     6     7     8     9
#  18105 13147 13096  9143  8936  6405  5414  4813   598   176


# 04. feature plots for annotating ----

# confirming clusters of cell types I originally subset 
astro_fts <- FeaturePlot(immune_vasc.obj, features = c("Aldh1l1", "Aqp4", "Gfap", "Slc7a10"))
micro_fts <- FeaturePlot(immune_vasc.obj, features = c("Tmem119", "P2ry12", "Aif1", "Trem2"))
BAM_fts <- FeaturePlot(immune_vasc.obj, features = c("Mrc1", "Mgl2"))
endo_fts <- FeaturePlot(immune_vasc.obj, features = c("Flt1", "Cldn5", "Emcn"))
fibro_fts <- FeaturePlot(immune_vasc.obj, features = c("Col1a2", "Lum"))
peri_fts <- FeaturePlot(immune_vasc.obj, features = c("Vtn", "Kcnj8", "Atp13a5"))
VLMC_fts <- FeaturePlot(immune_vasc.obj, features = c("Slc47a1", "Mgp", "Bgn"))
VSMC_fts <- FeaturePlot(immune_vasc.obj, features = c("Acta2", "Tagln"))
Tcell_fts <- FeaturePlot(immune_vasc.obj, features = c("Cd3d", "Gzmb", "Pdcd1"))

# other cell markers that should not be highly expressed in this obj
neu_fts <- FeaturePlot(immune_vasc.obj, features = c("Snap25", "Grin2a", "Grin2b"))
CPepend_fts <- FeaturePlot(immune_vasc.obj, features = c("Car12", "Ccdc153"))
oligoOPC_fts <- FeaturePlot(immune_vasc.obj, features = c("Ermn", "Cldn11", "Pdgfra"))


# save figures 
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_UMAPres04.png"), plot = UMAP_res04, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_imgdimplot_split_res04.png"), plot = imgdim_splitclusters_res04, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_imgdim_res04.png"), plot = imgdim_res04, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_ftplot.png"), plot = astro_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_micro_ftplot.png"), plot = micro_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_BAM_ftplot.png"), plot = BAM_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_oligoOPC_ftplot.png"), plot = oligoOPC_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_neu_ftplot.png"), plot = neu_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_CPepend_ftplot.png"), plot = CPepend_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_endo_ftplot.png"), plot = endo_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_fibro_ftplot.png"), plot = fibro_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_peri_ftplot.png"), plot = peri_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_VLMC_ftplot.png"), plot = VLMC_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_VSMC_ftplot.png"), plot = VSMC_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_Tcell_ftplot.png"), plot = Tcell_fts, path = figs_path)


# save fully processed obj with cluster info 
qs_save(immune_vasc.obj, paste0(output_path, date, "SpaNorm_RPCA_immunevasc_wclus.qs2"))


# read in obj
immune_vasc.obj <- qs_read(paste0(output_path, "20260814_SpaNorm_RPCA_immunevasc_wclus.qs2"))

Idents(immune_vasc.obj) <- "res04_clus"

# subset out clusters 8 bc it highly expresses many cell type  markers
immune_vasc.obj <- subset(immune_vasc.obj, idents = "8", invert = T)

dim(immune_vasc.obj)
#   480 79235

table(immune_vasc.obj$res04_clus, immune_vasc.obj$Age)
#    16   36   59   92
#0 4617 4632 4509 4347
#1 3559 3229 3155 3204
#2 3381 3370 3424 2921
#3 2144 2196 2196 2607
#4 2230 2049 2226 2431
#5 1449 1692 1622 1642
#6 1339 1410 1400 1265
#7 1207 1253 1191 1162
#8    0    0    0    0
#9   18   31   49   78

immune_vasc.obj <- FindSubCluster(immune_vasc.obj, cluster = "4", resolution = 0.2, graph.name = "integrated_snn")

Idents(immune_vasc.obj) <- "sub.cluster"

immune_vasc.obj$sub4_clus <- Idents(immune_vasc.obj)

# dimplot
UMAP_sub4 <- DimPlot(immune_vasc.obj, group.by = "sub4_clus", label = T, cols = "polychrome")

imgdim_sub4 <- ImageDimPlot(immune_vasc.obj, group.by = "sub4_clus", cols = "polychrome")

imgdim_split_sub4 <- ImageDimPlot(immune_vasc.obj, split.by = "sub4_clus", cols = "polychrome")

ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_sub4_UMAP.png"), plot = UMAP_sub4, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_sub4_imgdim.png"), plot = imgdim_sub4, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_sub4_imgdimsplit.png"), plot = imgdim_split_sub4, path = figs_path)


# find markers for sub 8 clusters
immunevascsub4.markers <- FindAllMarkers(immune_vasc.obj, group.by = "sub4_clus", only.pos = T, logfc.threshold = 0.25)

immunevascsub4.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_immunevascsub4 

write_csv(top5_immunevascsub4, paste0(cluster_markers, date, "SpaNorm_RPCA_immunevasc_sub4.csv"))

# save obj
qs_save(immune_vasc.obj, paste0(output_path, date, "SpaNorm_RPCA_immunevasc_sub4_clus.qs2"))

# read obj back in
immune_vasc.obj <- qs_read(paste0(output_path, "20260814_SpaNorm_RPCA_immunevasc_sub4_clus.qs2"))

# subcluster at lower resolution (res 0.2 for cluster 4 gave over 20 subclusters)
table(immune_vasc.obj$sub4_clus)
#  6     2     0     3     5     7   4_0     9   4_5   4_8   4_1     1   4_4
#5414 13096 18105  9143  6405  4813  2099   176   537     3  1973 13147   763
#4_2   4_3   4_6  4_16   4_9  4_10  4_11  4_17  4_19  4_12  4_13  4_14  4_20
#1749  1545   237     2     2     2     2     2     2     2     2     2     2
#4_21   4_7  4_15  4_18
#  2     4     2     2


Idents(immune_vasc.obj) <- "res04_clus"

immune_vasc.obj <- FindSubCluster(immune_vasc.obj, cluster = "4", resolution = 0.1, graph.name = "integrated_snn")

Idents(immune_vasc.obj) <- "sub.cluster"

immune_vasc.obj$sub4_res01_clus <- Idents(immune_vasc.obj)

# dimplot
UMAP_sub4 <- DimPlot(immune_vasc.obj, group.by = "sub4_res01_clus", label = T, cols = "polychrome")

imgdim_sub4 <- ImageDimPlot(immune_vasc.obj, group.by = "sub4_res01_clus", cols = "polychrome")

imgdim_split_sub4 <- ImageDimPlot(immune_vasc.obj, split.by = "sub4_res01_clus", cols = "polychrome")

ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_sub4_res01_UMAP.png"), plot = UMAP_sub4, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_sub4_res01_imgdim.png"), plot = imgdim_sub4, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_immunevasc_sub4_res01_imgdimsplit.png"), plot = imgdim_split_sub4, path = figs_path)


# find markers for sub 8 clusters
immunevascsub4.markers <- FindAllMarkers(immune_vasc.obj, group.by = "sub4_res01_clus", only.pos = T, logfc.threshold = 0.25)

immunevascsub4.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_immunevascsub4 

write_csv(top5_immunevascsub4, paste0(cluster_markers, date, "SpaNorm_RPCA_immunevasc_sub4_res01.csv"))

# save obj
qs_save(immune_vasc.obj, paste0(output_path, date, "SpaNorm_RPCA_immunevasc_sub4_clus.qs2"))

immune_vasc.obj <- qs_read(paste0(output_path, date, "SpaNorm_RPCA_immunevasc_sub4_clus.qs2"))

table(immune_vasc.obj$sub4_res01_clus)
#   6     2     0     3     5     7   4_1     9   4_3   4_5   4_2     1   4_0
# 5414 13096 18105  9143  6405  4813  2321   176  1293     3  1939 13147  3350
#  4_14   4_6   4_7   4_8  4_16   4_9  4_15  4_10  4_11  4_12  4_18   4_4  4_17
#    2     2     2     2     2     2     2     2     2     2     2     4     2
#4_13
#  2
 
Idents(immune_vasc.obj) <- "sub4_res01_clus"

# remove all sub4 res 01 clusters with small numbers of cells 
immune_vasc.obj <- subset(immune_vasc.obj, idents = c("0", "1", "2", "3", "4_0", "4_1", "4_2", "4_3", "5", "6", "7", "9"))

# rename idents and add to metadata ----
immune_vasc.obj <- RenameIdents(immune_vasc.obj, 
             "0" = "Endothelial cells", 
             "1" = "Astrocytes", 
             "2" = "Astrocytes", 
             "3" = "Microglia", 
             "4_0" = "Astrocytes",
             "4_1" = "Endothelial cells", 
             "4_2" = "VLMC",
             "4_3" = "BAM", 
             "5" = "Astrocytes", 
             "6" = "VSMC", 
             "7" = "Pericytes", 
             "9" = "T cells"
             )


immune_vasc.obj$immune_vasc_ann <- Idents(immune_vasc.obj)

# save object annotated obj ----
qs_save(immune_vasc.obj, paste0(output_path, date, "SpaNorm_RPCA_immunevasc_annotated.qs2"))
