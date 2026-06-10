########################################################################
# Name: CCL_ST_WTAgingXMet_SCT_03.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 03. 
#             - 01
#             - 02 
#             - 03
#             - 04  
# Input Files:  
# Output Files: 
# Date created: 4/15/26
# Last updated: 4/24/26
# Author: Chloe Lucido
########################################################################

library(spacexr)
library(Seurat)
library(SeuratDisk)
library(future)
library(ggplot2)
library(arrow)
library(hdf5r)
library(presto)
library(glmGamPoi)
library(readr)
library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)
library(patchwork)
library(sceasy)
library(reticulate)
library(SPLIT)
library(qs2) 
library(RColorBrewer)
library(Polychrome)
library(purrr)
library(readxl)
library(BPCells)
options(future.globals.maxSize = 100 *1024^3)

# set seed
set.seed(42)


WTmerged.obj <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260421_mergedobj_SCT.qs2")

UMAP_v1 <- DimPlot(WTmerged.obj, label = T, cols = "polychrome")

# findallmarkers 

WTmerged.obj.markers <- FindAllMarkers(WTmerged.obj, group.by = "seurat_clusters")

WTmerged.obj.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_WT_xenium

write_csv(top5_WT_xenium, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/cluster_markers/20260422_SCT_top5_WT_xenium_res05.csv")

table(WTmerged.obj$seurat_clusters)

ggsave(filename = "20260422_allcells_SCT_UMAP_res05.png", plot = UMAP_v1, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs/", width = 9, height = 10)

# saving res0.5 clusters in metadata 
WTmerged.obj$res05_clusters <- Idents(WTmerged.obj)

# changing res to 0.4 
WTmerged.obj <- FindClusters(WTmerged.obj, resolution = 0.4, graph.name = "SCT_snn")

UMAP_v2 <- DimPlot(WTmerged.obj, label = T, cols = "polychrome")

WTmerged.obj.markers <- FindAllMarkers(WTmerged.obj, group.by = "seurat_clusters")

WTmerged.obj.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_WT_xenium04

write_csv(top5_WT_xenium04, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/cluster_markers/20260422_SCT_top5_WT_xenium_res04.csv")

ggsave(filename = "20260422_allcells_SCT_UMAP_res04.png", plot = UMAP_v2, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs/", width = 9, height = 10)

# saving res 0.4 clusters in metadata 
WTmerged.obj$res04_clusters <- Idents(WTmerged.obj)

# changing res to 0.3
WTmerged.obj <- FindClusters(WTmerged.obj, resolution = 0.3, graph.name = "SCT_snn")

UMAP_v3 <- DimPlot(WTmerged.obj, label = T, cols = "polychrome")

WTmerged.obj.markers <- FindAllMarkers(WTmerged.obj, group.by = "seurat_clusters")

WTmerged.obj.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_WT_xenium03

write_csv(top5_WT_xenium03, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/cluster_markers/20260422_SCT_top5_WT_xenium_res03.csv")

ggsave(filename = "20260422_allcells_SCT_UMAP_res03.png", plot = UMAP_v3, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs/", width = 9, height = 10)

# saving res 0.3 clusters in metadata 
WTmerged.obj$res03_clusters <- Idents(WTmerged.obj)

# save SCT obj w/ clusters metadata 
qs_save(WTmerged.obj, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260422_mergedobj_SCT_clusteredres03-05.qs2")


# annotate clusters 

FeaturePlot(WTmerged.obj, features = c("Aldh1l1", "Gfap", "Slc7a10", "Aqp4")) # astro

FeaturePlot(WTmerged.obj, features = c("Tmem119", "Aif1", "Trem2", "P2ry12")) # microglia

FeaturePlot(WTmerged.obj, features = c("Cst7", "Clec7a"))

FeaturePlot(WTmerged.obj, features = c("Ermn", "Cldn11", "Mog")) # oligo

FeaturePlot(WTmerged.obj, features = c("Ccdc153", "Dnah11", "Tmem212")) # epend

FeaturePlot(WTmerged.obj, features = c("Pdgfra", "Tnr")) # OPC

FeaturePlot(WTmerged.obj, features = c("Mgp", "Bgn", "Slc47a1")) # VLMC

FeaturePlot(WTmerged.obj, features = c("Acta2", "Tagln")) # VSMC

FeaturePlot(WTmerged.obj, features = c("Vtn", "Kcnj8", "Atp13a5")) # pericytes

FeaturePlot(WTmerged.obj, features = c("Flt1", "Emcn", "Cldn5")) # endo

FeaturePlot(WTmerged.obj, features = c("Car12", "Ttr")) # CP

FeaturePlot(WTmerged.obj, features = c("Col1a2", "Lum")) # fibroblast

FeaturePlot(WTmerged.obj, features = c("Dcx", "Pax6")) # NPC

FeaturePlot(WTmerged.obj, features = c("Mgl2", "Mrc1")) # BAM

FeaturePlot(WTmerged.obj, features = c("Cd3d", "Gzmb", "Pdcd1")) # T cell

FeaturePlot(WTmerged.obj, features = c("Snap25", "Grin2a")) # pan-neu

FeaturePlot(WTmerged.obj, features = c("Otof", "Lypd1", "Pou3f1", "Syt6")) # excitatory neu

FeaturePlot(WTmerged.obj, features = c("Lamp5", "Pvalb", "Calb1", "Sst"))




