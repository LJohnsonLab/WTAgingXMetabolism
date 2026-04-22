########################################################################
# Name: CCL_ST_WTAgingXMet_02.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 02. merging all objs and processing (NormalizeData->FindVariableFeatures->ScaleData) 
#             - 01 
#             - 02 
#             - 03  
# Input Files:  
# Final Output Files: 
# Date created: 4/15/26
# Last updated: 4/22/26
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

# 01. read in individual objs ----

xenium.objWT_999 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0021999.qs2")
xenium.objWT_573 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0022573.qs2")
xenium.objWT_991 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260420_WT_slide0021991.qs2")
xenium.objWT_998 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0021998.qs2")
xenium.objWT_039 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0069039.qs2")
xenium.objWT_045 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0069045.qs2")
xenium.objWT_118 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260420_WT_slide0069118.qs2")


# 02. Merge objects and saving merged file ----

## 02a. make list of all individual WT subsetted objs ----
xenium_list <- list(xenium.objWT_999, 
                 xenium.objWT_573, 
                 xenium.objWT_991, 
                 xenium.objWT_998, 
                 xenium.objWT_039, 
                 xenium.objWT_045, 
                 xenium.objWT_118)

## 02b. merge objects using the list created above ----
WTmerged.obj <- merge(xenium_list[[1]], y = xenium_list[-1])

# SANITY CHECK
dim(WTmerged.obj)
# 480 1619235

## 02c. Join layers ----
WTmerged.obj <- JoinLayers(WTmerged.obj, add.prefix = FALSE)

## 02d. subset nCount > 5 ----
WTmerged.obj <- subset(WTmerged.obj, subset = nCount_Xenium > 5)

# SANITY CHECK
dim(WTmerged.obj)
# 480 1609408

# 03. save merged object ----

qs_save(WTmerged.obj, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260420_mergedobj_no_metadata.qs2")

# 04. Attaching metadata from xlsx file ----

## 04a. load merged object ----

WTmerged.obj <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260420_mergedobj_no_metadata.qs2")

## 04b. read in metadata excel sheet ----
WT_metadata <- data.frame(
  read_excel("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/metadata/WT_AgingXMetabolism_Xenium_metadata.xlsx"))

## 04c. extract metadata from merged seurat obj ----
meta <- WTmerged.obj@meta.data

## 04d. convert rownames to cell id ----
meta <- meta %>%
  rownames_to_column(var = "cell_id")

## 04e. join metadata excel sheet to merged obj metadata 
meta <- meta %>%
  left_join(WT_metadata, by = "sample_ID"
            #, keep = T
  )

## 04f. revert cell id back to columns ----
meta <- column_to_rownames(meta, var = "cell_id")

## 04g. attach metadata back to merged obj ----
WTmerged.obj@meta.data <- meta


# SANITY CHECK
dim(WTmerged.obj)

table(WTmerged.obj$sample_ID, WTmerged.obj$slide)

table(WTmerged.obj$Age, WTmerged.obj$Sex)

table(WTmerged.obj$sample_ID, WTmerged.obj$Age)

# 05. save merged obj ----
qs_save(WTmerged.obj, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260420_mergedobjv1_01.qs2")



# 06. read in merged obj ----
WT_xenium <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260420_mergedobjv1_01.qs2")


# 07. Update merged obj ----

WT_xenium <- UpdateSeuratObject(WT_xenium)


VlnPlot(WT_xenium, features = c("nFeature_Xenium", "nCount_Xenium"), 
                  raster = F,  
                  ncol = 2, pt.size = 0)




# 08. Standard Seurat pipeline ----

# 4/17/26: everytime I try to do SCTransform it won't work bc of memory usage

## 08a. Normalize ----
WT_xenium <- NormalizeData(WT_xenium)

## 08b. FindVariableFeatures ----
WT_xenium <- FindVariableFeatures(WT_xenium)

## 08c. ScaleData ----
WT_xenium <- ScaleData(WT_xenium, vars.to.regress = "slide_ID")

## 08d. RunPCA ----
WT_xenium <- RunPCA(WT_xenium, npcs = 50, features = rownames(WT_xenium))

PCAPlot(WT_xenium, group.by = "Sex")

# elbow plot to determine number of PCs to use 
ElbowPlot(WT_xenium, ndims = 50, reduction = "pca")

## 08e. RunUMAP ----
WT_xenium <- RunUMAP(WT_xenium, dims = 1:30)

## 08f. FindNeighbors 
WT_xenium <- FindNeighbors(WT_xenium, reduction = "pca", dims = 1:30)

## 08g. FindClusters 
WT_xenium <- FindClusters(WT_xenium, resolution = 0.3, graph.name = "Xenium_snn")

## 08h. Visualize UMAP with clusters 
DimPlot(WT_xenium, label = T, cols = "polychrome")

# UMAP SAVED AS 20260420_allcells_UMAP_normalize_scale_res03.png in /figs/UMAPs


## 08i. FindAllMarkers 
WT_xenium.markers <- FindAllMarkers(WT_xenium, group.by = "seurat_clusters")

WT_xenium.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_WT_xenium

write_csv(top5_WT_xenium, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/cluster_markers/20260420_top5_WT_xenium_res03.csv")

# annotate clusters 

FeaturePlot(WT_xenium, features = c("Aldh1l1", "Gfap", "Slc7a10", "Aqp4")) # astro

FeaturePlot(WT_xenium, features = c("Tmem119", "Aif1", "Trem2", "P2ry12")) # mircoglia

FeaturePlot(WT_xenium, features = c("Cst7", "Clec7a"))

FeaturePlot(WT_xenium, features = c("Ermn", "Cldn11", "Mog")) # oligo

FeaturePlot(WT_xenium, features = c("Ccdc153", "Dnah11", "Tmem212")) # epend

FeaturePlot(WT_xenium, features = c("Pdgfra", "Tnr")) # OPC

FeaturePlot(WT_xenium, features = c("Mgp", "Bgn", "Slc47a1")) # VLMC

FeaturePlot(WT_xenium, features = c("Acta2", "Tagln")) # VSMC

FeaturePlot(WT_xenium, features = c("Vtn", "Kcnj8", "Atp13a5")) # pericytes

FeaturePlot(WT_xenium, features = c("Flt1", "Emcn", "Cldn5")) # endo

FeaturePlot(WT_xenium, features = c("Car12", "Ttr")) # CP

FeaturePlot(WT_xenium, features = c("Col1a2", "Lum")) # fibroblast

FeaturePlot(WT_xenium, features = c("Dcx", "Pax6")) # NPC

FeaturePlot(WT_xenium, features = c("Mgl2", "Mrc1")) # BAM

FeaturePlot(WT_xenium, features = c("Cd3d", "Gzmb", "Pdcd1")) # T cell

FeaturePlot(WT_xenium, features = c("Snap25", "Grin2a")) # pan-neu

FeaturePlot(WT_xenium, features = c("Otof", "Lypd1", "Pou3f1", "Syt6")) # excitatory neu

FeaturePlot(WT_xenium, features = c("Lamp5", "Pvalb", "Calb1", "Sst"))

# 09. UPDATE THIS ----
qs_save(WT_xenium, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Files/") # UPDATE THIS
