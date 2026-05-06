########################################################################
# Name: CCL_ST_WTAgingXMet_04.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 04. Subset, subcluster and annotate neuronal subtypes  
#             - 01 read in annotated obj from part 3
#             - 02 subset neurons
#             - 03 re-process neurons-only obj
#             - 04 plot neuronal subtype markers for annotation (feature plots)
#             - 05 save obj before removing genes => 20260430_neu_res03.qs2
#             - 06 remove genes 
#             - 07 re-process obj 
#             - 08 save obj with genes removed 
# Input Files: 20260428_mergedobj_gen_annotations_03.qs2  
# Final Output Files: 20260505_neu_res03_genesremoved.qs2
# Date created: 4/28/26
# Last updated: 5/5/26
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


# 02. subset neurons from full obj ----
neurons.obj <- subset(WTmerged.obj, idents = "Neuron") 

# SANITY CHECK
DimPlot(neurons.obj)

# 03. re-process neuron obj ----

## 03a. Normalize ----
neurons.obj <- NormalizeData(neurons.obj)

## 03b. FindVariableFeatures ----
neurons.obj <- FindVariableFeatures(neurons.obj)

## 03c. ScaleData ----
neurons.obj <- ScaleData(neurons.obj)

## 03d. RunPCA ----
neurons.obj <- RunPCA(neurons.obj, npcs = 40, features = rownames(neurons.obj))

PCAPlot(neurons.obj, group.by = "Sex")

# elbow plot to determine number of PCs to use 
ElbowPlot(neurons.obj, ndims = 40, reduction = "pca")

## 03e. RunUMAP ----
neurons.obj <- RunUMAP(neurons.obj, dims = 1:25)

## 03f. FindNeighbors ----
neurons.obj <- FindNeighbors(neurons.obj, reduction = "pca", dims = 1:25)

## 03g. FindClusters ----
neurons.obj <- FindClusters(neurons.obj, resolution = 0.3, graph.name = "Xenium_snn")

## 03h. Visualize UMAP with clusters ----
UMAP_res03 <- DimPlot(neurons.obj, label = T, cols = "polychrome")

imgdim_clusters_res03 <- ImageDimPlot(neurons.obj, split.by = "seurat_clusters", fov = "fov")

imgdim_res03 <- ImageDimPlot(neurons.obj, group.by = "seurat_clusters", cols = "polychrome")


## 03i. find markers ----
neurons.markers <- FindAllMarkers(neurons.obj, group.by = "seurat_clusters", only.pos = T, logfc.threshold = 0.25)

neurons.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_neurons

write_csv(top5_neurons, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/cluster_markers/20260428_top5_res03_neu.csv")

ggsave(filename = "20260429_neu_UMAP_res03.png", plot = UMAP_res03, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs/", width = 11, height = 10)
ggsave(filename = "20260429_neu_imgdim_clusters_res03.png", plot = imgdim_clusters_res03, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/ImageDimPlot/", width = 11, height = 10)
ggsave(filename = "20260429_neu_imgdim_res03.png", plot = imgdim_res03, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/ImageDimPlot/", width = 11, height = 10)

# save clusters in metadata 
neurons.obj$neu_res03_clusters <- Idents(neurons.obj)

# 04. feature plots ----

FeaturePlot(neurons.obj, features = c("Snap25", "Grin2a", "Grin2b")) # pan-neu

FeaturePlot(neurons.obj, features = c("Fezf2", "Slc17a7", "Slc17a6", "Htr2c")) # general excitatory neu

FeaturePlot(neurons.obj, features = c("Cux2", "Otof", "Stard8", "Lypd1")) # L2/3 IT (exc) neu

FeaturePlot(neurons.obj, features = "Rapgef3") # L5 PT (exc) neu

FeaturePlot(neurons.obj, features = c("Rorb", "Tcap", "Rspo1", "Whrn")) # L5 IT (exc) neu

FeaturePlot(neurons.obj, features = c("Tunar", "Osr1", "Oprk1")) # L6 IT (exc) neu

FeaturePlot(neurons.obj, features = "Pou3f1") # L5 ET (exc) neu

FeaturePlot(neurons.obj, features = c("Tshz2", "Tox2")) # L5/6 NP (exc) neu

FeaturePlot(neurons.obj, features = c("Syt6", "Trh")) # L6 CT (exc) neu

FeaturePlot(neurons.obj, features = c("Gad1", "Lhx6", "Cnr1", "Serpinf1", "Slc32a1")) # general inhibitory neu

FeaturePlot(neurons.obj, features = c("Lamp5", "Ndnf")) # LAMP5 (inh) neu

FeaturePlot(neurons.obj, features = c("Pvalb", "Nts", "Tac1")) # PVALB (inh) neu

FeaturePlot(neurons.obj, features = "Sncg") # SNCG (inh) neu

FeaturePlot(neurons.obj, features = c("Sst", "Calb1")) # SST (inh) neu

FeaturePlot(neurons.obj, features = c("Vip", "Calb2", "Pthlh", "Crh")) # VIP (inh) neu

# 05. save obj before removing genes ----
qs_save(neurons.obj, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/Neuron_objs/20260430_neu_res03.qs2")

# LOAD IN OBJ ----

neurons.obj <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/Neuron_objs/20260430_neu_res03.qs2")

# 06. remove genes ----
# removing cell-type markers genes and other immune genes that should not be expressed by neurons 

## 06a. read in xlsx with genes to be removed ----
genes_to_remove <- data.frame(
  read_excel("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/genes_removed/Genes_removed_from_neurons.xlsx"))


# genes are stored as genes_to_remove$Genes

# SANITY CHECK: dimensions before removing genes 
dim(neurons.obj)
# 480 735213

## 06b. subset obj ----

# create list of genes to keep 
genes_to_keep <- setdiff(rownames(neurons.obj), genes_to_remove$Genes)

# subset obj to only contain genes to keep
neurons.obj_v1 <- subset(neurons.obj, features = genes_to_keep)

# SANITY CHECK 
dim(neurons.obj_v1)
# 291 735213

# 07. re-process ----

## 07a. Normalize ----
neurons.obj_v1 <- NormalizeData(neurons.obj_v1)

## 07b. FindVariableFeatures ----
neurons.obj_v1 <- FindVariableFeatures(neurons.obj_v1)

## 07c. ScaleData ----
neurons.obj_v1 <- ScaleData(neurons.obj_v1)

## 07d. RunPCA ----
neurons.obj_v1 <- RunPCA(neurons.obj_v1, npcs = 30, features = rownames(neurons.obj_v1))

PCAPlot(neurons.obj_v1, group.by = "Sex")

# elbow plot to determine number of PCs to use 
ElbowPlot(neurons.obj_v1, ndims = 30, reduction = "pca")

## 07e. RunUMAP ----
neurons.obj_v1 <- RunUMAP(neurons.obj_v1, dims = 1:25)

## 07f. FindNeighbors ----
neurons.obj_v1 <- FindNeighbors(neurons.obj_v1, reduction = "pca", dims = 1:25)

## 07g. FindClusters ----
neurons.obj_v1 <- FindClusters(neurons.obj_v1, resolution = 0.3, graph.name = "Xenium_snn")


UMAP_res03_v1 <- DimPlot(neurons.obj_v1, group.by = "seurat_clusters", cols = "polychrome")

imgdim_clusters_res03_v1 <- ImageDimPlot(neurons.obj_v1, split.by = "seurat_clusters", cols = "polychrome")
imgdim_res03_v1 <- ImageDimPlot(neurons.obj_v1, group.by = "seurat_clusters", cols = "polychrome")

ggsave(filename = "20260505_neu_genesremovedUMAP_res03.png", plot = UMAP_res03_v1, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs/", width = 11, height = 10)
ggsave(filename = "20260505_neu_genesremovedimgdim_clusters_res03.png", plot = imgdim_clusters_res03_v1, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/ImageDimPlot/", width = 11, height = 10)
ggsave(filename = "20260505_neu_genesremovedimgdim_res03.png", plot = imgdim_res03_v1, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/ImageDimPlot/", width = 11, height = 10)


neurons.markers_v1 <- FindAllMarkers(neurons.obj_v1, group.by = "seurat_clusters", only.pos = T, logfc.threshold = 0.25)

neurons.markers_v1 |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_neurons_v1

write_csv(top5_neurons_v1, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/cluster_markers/20260505_top5_res03_genesremovedneu.csv")

table(neurons.obj_v1$seurat_clusters)

FeaturePlot(neurons.obj_v1, features = c("Snap25", "Grin2a", "Grin2b")) # pan-neu

FeaturePlot(neurons.obj_v1, features = c("Fezf2", "Slc17a7", "Slc17a6", "Htr2c")) # general excitatory neu

FeaturePlot(neurons.obj_v1, features = c("Cux2", "Otof", "Stard8", "Lypd1")) # L2/3 IT (exc) neu

FeaturePlot(neurons.obj_v1, features = "Rapgef3") # L5 PT (exc) neu

FeaturePlot(neurons.obj_v1, features = c("Rorb", "Tcap", "Rspo1", "Whrn")) # L5 IT (exc) neu

FeaturePlot(neurons.obj_v1, features = c("Tunar", "Osr1", "Oprk1")) # L6 IT (exc) neu

FeaturePlot(neurons.obj_v1, features = "Pou3f1") # L5 ET (exc) neu

FeaturePlot(neurons.obj_v1, features = c("Tshz2", "Tox2")) # L5/6 NP (exc) neu

FeaturePlot(neurons.obj_v1, features = c("Syt6", "Trh")) # L6 CT (exc) neu

FeaturePlot(neurons.obj_v1, features = c("Gad1", "Lhx6", "Cnr1", "Serpinf1", "Slc32a1")) # general inhibitory neu

FeaturePlot(neurons.obj_v1, features = c("Lamp5", "Ndnf")) # LAMP5 (inh) neu

FeaturePlot(neurons.obj_v1, features = c("Pvalb", "Nts", "Tac1")) # PVALB (inh) neu

FeaturePlot(neurons.obj_v1, features = "Sncg") # SNCG (inh) neu

FeaturePlot(neurons.obj_v1, features = c("Sst", "Calb1")) # SST (inh) neu

FeaturePlot(neurons.obj_v1, features = c("Vip", "Calb2", "Pthlh", "Crh")) # VIP (inh) neu

## 7h. save clusters to metadata ----

Idents(neurons.obj_v1)

neurons.obj_v1$neu_res03_clusters_v1 <- Idents(neurons.obj_v1)


# 08. save object w/ genes removed ----

qs_save(neurons.obj_v1, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/Neuron_objs/20260505_neu_res03_genesremoved.qs2")
