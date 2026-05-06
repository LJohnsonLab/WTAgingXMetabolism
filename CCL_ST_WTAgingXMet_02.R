########################################################################
# Name: CCL_ST_WTAgingXMet_02.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 02. merging all objs and processing (NormalizeData->FindVariableFeatures->ScaleData) 
#             - 01 read in individual files 
#             - 02 merge objects
#             - 03 save merged obj => 20260420_mergedobj_no_metadata.qs2
#             - 04 attach metadata (from xlsx) to merged obj 
#             - 05 save merged obj with metadata => 20260420_mergedobjv1_01.qs2
#             - 06 update merged obj
#             - 07 add BPCells backed matrix
#             - 08 standard seurat pipeline (normalize and scale instead of SCTransform)
# Input Files: individual obj files from part 1  
# Final Output Files: 
# Date created: 4/15/26
# Last updated: 4/24/26
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



# 06. Update merged obj ----

## 06a. read in merged obj ----
WTmerged.obj <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260420_mergedobjv1_01.qs2")


## 06b. Update merged obj ----

WTmerged.obj <- UpdateSeuratObject(WTmerged.obj)


VlnPlot(WTmerged.obj, features = c("nFeature_Xenium", "nCount_Xenium"), 
                  raster = F,  
                  ncol = 2, pt.size = 0)




# 07. Add BPCells backed matrix to reduce memory ----

## 07a. create BPCells directory for raw counts to be stored in ----
write_matrix_dir(
  mat = convert_matrix_type( # ensure the counts are stored as integers and not doubles (aka decimals) to compress efficiently
    LayerData(WTmerged.obj, layer = "counts"), # RAW COUNTS (not normalized or scaled)
    type = "uint32_t"), # unsigned 32-bit integer
  dir = "./bpcells/mergedobj_counts"
)


## 07b. open BPCells directory ----
bp_counts <- open_matrix_dir("./bpcells/mergedobj_counts")

## 07c. change dgcMatrix to BPCells matrix ----
WTmerged.obj[["Xenium"]] <- CreateAssay5Object(counts = bp_counts)


# 08. Standard Seurat pipeline ----

# 4/17/26: everytime I try to do SCTransform it won't work bc of memory usage

## 08a. Normalize ----
WTmerged.obj <- NormalizeData(WTmerged.obj)

## 08b. FindVariableFeatures ----
WTmerged.obj <- FindVariableFeatures(WTmerged.obj)

## 08c. ScaleData ----
WTmerged.obj <- ScaleData(WTmerged.obj, vars.to.regress = "slide_ID") # regressed out slide_ID to prevent batch effects

## 08d. RunPCA ----
WTmerged.obj <- RunPCA(WTmerged.obj, npcs = 50, features = rownames(WTmerged.obj))

PCAPlot(WTmerged.obj, group.by = "Sex")

# elbow plot to determine number of PCs to use 
ElbowPlot(WTmerged.obj, ndims = 50, reduction = "pca")

## 08e. RunUMAP ----
WTmerged.obj <- RunUMAP(WTmerged.obj, dims = 1:30)

## 08f. FindNeighbors ----
WTmerged.obj <- FindNeighbors(WTmerged.obj, reduction = "pca", dims = 1:30)

## 08g. FindClusters ----
WTmerged.obj <- FindClusters(WTmerged.obj, resolution = 0.4, graph.name = "Xenium_snn")

## 08h. Visualize UMAP with clusters ----
DimPlot(WTmerged.obj, label = T, cols = "polychrome")


# 09. Save normalized and scaled obj ----
qs_save(WTmerged.obj, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260424_mergedobj_res04_02.qs2") 
