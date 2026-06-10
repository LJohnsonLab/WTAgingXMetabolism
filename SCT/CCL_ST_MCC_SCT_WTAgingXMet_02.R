########################################################################
# Name: CCL_ST_MCC_SCT_WTAgingXMet_02.R **ONLY TO BE RUN ON SUPERCOMPUTER (MCC)**
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - MCC 02. merging objs, attaching metadata, and running SCTransform using MCC 
#             - 01 merge objs and subset nCount>5
#             - 02 save merged obj (no metadata) => 20260421_mergedobj_no_metadata.qs2
#             - 03 attach metadata (xlsx) to merged obj
#             - 04 Processing (SCTransform)
#             - 05 save SCTransformed merged obj => 20260421_mergedobj_SCT.qs2
# Input Files: individual objs from part 1
# Output Files: SCTransformed merged obj => 20260421_mergedobj_SCT.qs2
# Date created: 4/15/26
# Last updated: 6/10/26 (updated file name in header)
# Author: Chloe Lucido
########################################################################


# Auto-install missing packages
required_packages <- c("Seurat", "future", "ggplot2", "readr", "dplyr", 
                       "data.table","tidyr", "tibble", "stringr", "forcats", "lubridate",  "readxl", "patchwork", 
                       "qs2", "RColorBrewer", "Polychrome", "purrr",
                       "presto", "glmGamPoi")

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

# ONLY TO BE RUN ON MCC 

# Use the cores you requested
plan("multicore", workers = 16)
options(future.globals.maxSize = 240 * 1024^3)  # slightly under your --mem

# set seed
set.seed(42)

# to prevent this error: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so.
options(future.rng.onMisuse = "resolve")
options(future.seed = TRUE)


# 01. Merge objects and saving merged file ----

## 01a. read in objs ----

xenium.objWT_999 <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/20260416_WT_slide0021999.qs2")
xenium.objWT_573 <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/20260416_WT_slide0022573.qs2")
xenium.objWT_991 <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/20260420_WT_slide0021991.qs2")
xenium.objWT_998 <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/20260416_WT_slide0021998.qs2")
xenium.objWT_039 <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/20260416_WT_slide0069039.qs2")
xenium.objWT_045 <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/20260416_WT_slide0069045.qs2")
xenium.objWT_118 <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/20260420_WT_slide0069118.qs2")



## 01b. make list of all individual WT subsetted objs ----
xenium_list <- list(xenium.objWT_999, 
                 xenium.objWT_573, 
                 xenium.objWT_991, 
                 xenium.objWT_998, 
                 xenium.objWT_039, 
                 xenium.objWT_045, 
                 xenium.objWT_118)

## 01c. merge objets using the list created above ----
WTmerged.obj <- merge(xenium_list[[1]], y = xenium_list[-1])

# SANITY CHECK
dim(WTmerged.obj)
# 480 1619235

## 01d. Join layers ----
WTmerged.obj <- JoinLayers(WTmerged.obj, add.prefix = FALSE)

## 01e. subset nCount > 5 ----
WTmerged.obj <- subset(WTmerged.obj, subset = nCount_Xenium > 5)

# SANITY CHECK
dim(WTmerged.obj)
# 480 1609408

# 02. save merged object ----

qs_save(WTmerged.obj, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260421_mergedobj_no_metadata.qs2") 



# 03. Attaching metadata from xlsx file ----


## 03a. read in metadata excel sheet ----
WT_metadata <- data.frame(
  read_excel("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/WT_AgingXMetabolism_Xenium_metadata.xlsx")) 

## 03b. extract metadata from merged seurat obj ----
meta <- WTmerged.obj@meta.data

## 03c. convert rownames to cell id ----
meta <- meta %>%
  rownames_to_column(var = "cell_id")

## 03d. join metadata excel sheet to merged obj metadata ----
meta <- meta %>%
  left_join(WT_metadata, by = "sample_ID"
            #, keep = T
  )

## 03e. revert cell id back to columns ----
meta <- column_to_rownames(meta, var = "cell_id")

## 3f. attach metadata back to merged obj ----
WTmerged.obj@meta.data <- meta


# 04. Processing ----

WTmerged.obj <- UpdateSeuratObject(WTmerged.obj)


## 04a. Standard Seurat pipeline ----


WTmerged.obj <- SCTransform(WTmerged.obj, assay = "Xenium")

WTmerged.obj <- RunPCA(WTmerged.obj, npcs = 50, features = rownames(WTmerged.obj))

WTmerged.obj <- RunUMAP(WTmerged.obj, dims = 1:30)

WTmerged.obj <- FindNeighbors(WTmerged.obj, reduction = "pca", dims = 1:30)

WTmerged.obj <- FindClusters(WTmerged.obj, resolution = 0.5, graph.name = "SCT_snn")


# 05. Save SCTransformed merged obj ----
qs_save(WTmerged.obj, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260421_mergedobj_SCT.qs2") 

