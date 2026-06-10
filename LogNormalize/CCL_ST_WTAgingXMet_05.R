########################################################################
# Name: CCL_ST_WTAgingXMet_05.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 04. Subset, subcluster and annotate vascular mural cell types  
#             - 01 read in annotated obj from part 3
#             - 02 subset vascular mural population
#             - 03 re-process vascular mural-only obj
#             - 04 plot vascular mural subtype markers for annotation (feature plots)
#             - 05 save obj before removing genes => 
#             - 06 remove genes 
#             - 07 re-process obj 
#             - 08 save obj with genes removed 
# Input Files: 20260428_mergedobj_gen_annotations_03.qs2  
# Final Output Files: 
# Date created: 5/5/26
# Last updated: 6/10/26
# Author: Chloe Lucido
########################################################################

# Load Libraries ----
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

# 01. read in final annotated obj from part 3 ----
WTmerged.obj <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260428_mergedobj_gen_annotations_03.qs2")


# 02. subset vascular mural from full obj ----
mural.obj <- subset(WTmerged.obj, idents = "Vascular Mural") 

# SANITY CHECK
DimPlot(mural.obj)

# SANITY CHECK
dim(mural.obj)
# 480 99376

# 03. re-process neuron obj ----

## 03a. Normalize ----
mural.obj <- NormalizeData(mural.obj)

## 03b. FindVariableFeatures ----
mural.obj <- FindVariableFeatures(mural.obj)

## 03c. ScaleData ----
mural.obj <- ScaleData(mural.obj)

## 03d. RunPCA ----
mural.obj <- RunPCA(mural.obj, npcs = 30, features = rownames(mural.obj)) 

PCAPlot(mural.obj, group.by = "Age")

# elbow plot to determine number of PCs to use 
ElbowPlot(mural.obj, ndims = 30, reduction = "pca")

## 03e. RunUMAP ----
mural.obj <- RunUMAP(mural.obj, dims = 1:20)

## 03f. FindNeighbors ----
mural.obj <- FindNeighbors(mural.obj, reduction = "pca", dims = 1:20)

## 03g. FindClusters ----
mural.obj <- FindClusters(mural.obj, resolution = 0.2, graph.name = "Xenium_snn")

## 03h. Visualize UMAP with clusters ----
mural_UMAP_res02 <- DimPlot(mural.obj, label = T, cols = "polychrome") 

mural_imgdim_clusters_res02 <- ImageDimPlot(mural.obj, split.by = "seurat_clusters", fov = "fov") 

mural_imgdim_res02 <- ImageDimPlot(mural.obj, group.by = "seurat_clusters", cols = "polychrome") 

ggsave(filename = "20260514_mural_UMAP_res02.png", plot = mural_UMAP_res02, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs/", width = 11, height = 10) 
ggsave(filename = "20260514_mural_imgdim_clusters_res02.png", plot = mural_imgdim_clusters_res02, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/ImageDimPlot/", width = 11, height = 10) 
ggsave(filename = "20260514_mural_imgdim_res02.png", plot = mural_imgdim_res02, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/ImageDimPlot/", width = 11, height = 10) 


table(mural.obj$seurat_clusters)

## 03i. find markers ----
mural.markers <- FindAllMarkers(mural.obj, group.by = "seurat_clusters", only.pos = T, logfc.threshold = 0.25)

mural.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  ungroup() -> top10_mural

write_csv(top5_mural, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/cluster_markers/20260506_top5_res02_mural.csv") 
write_csv(top10_mural, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/cluster_markers/20260514_top10_res02_mural.csv")

FeaturePlot(mural.obj, features = c("Mgp", "Bgn", "Slc47a1")) # VLMC
FeaturePlot(mural.obj, features = c("Acta2", "Tagln")) # VSMC
FeaturePlot(mural.obj, features = c("Mgl2", "Mrc1")) # BAM

# checking for contamination
FeaturePlot(mural.obj, features = c("Aqp4", "Gfap", "Aldoc", "Slc7a10")) # astro markers
FeaturePlot(mural.obj, features = c("Snap25", "Grin2a", "Cnr1", "Otof")) # various neuronal markers 
FeaturePlot(mural.obj, features = c("Mog", "Cldn11", "Ermn")) # oligos
FeaturePlot(mural.obj, features = c("Vtn", "Kcnj8")) # pericytes
FeaturePlot(mural.obj, features = c("Emcn", "Cldn5", "Flt1")) # endothelial cells



# save clusters in metadata 
#mural.obj$_ <- Idents(mural.obj) # CHANGE THIS

# 04. save obj before removing genes ----
qs_save(mural.obj, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/vascmural_objs/20260513_mural_res02.qs2")

