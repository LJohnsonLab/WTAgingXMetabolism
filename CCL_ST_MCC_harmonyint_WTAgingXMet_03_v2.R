########################################################################
# Name: CCL_ST_MCC_harmonyint_WTAgingXMet_03_v2.R **ONLY TO BE RUN ON SUPERCOMPUTER (MCC)**
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
# Date created: 5/13/26
# Last updated: 5/26/26
# Author: Chloe Lucido
# Version history:
#       -v2 following seurat v5 integrative analysis vignette https://satijalab.org/seurat/articles/seurat5_integration 
########################################################################


# Auto-install missing packages
required_packages <- c("Seurat", "future", "ggplot2", "readr", "dplyr",
                        "data.table","tidyr", "tibble", "stringr", "forcats", "lubridate",  "readxl", "patchwork",
                        "qs2", "RColorBrewer", "Polychrome", "purrr",
                        "presto", "glmGamPoi", "harmony")

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
library(harmony)
# ONLY TO BE RUN ON MCC

# Use the cores you requested
plan("multicore", workers = 8)
options(future.globals.maxSize = 650 * 1024^3)  # slightly under your --mem


# set seed
set.seed(42)

# to prevent this error: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so.
options(future.rng.onMisuse = "resolve")
options(future.seed = TRUE)
 

 
# 01. read in merged obj (not processed) no  metadata attached ----
# 5/15 objs: accidentally used 20260421_mergedobj_SCT.qs2 instead of 20260421_mergedobj_no_metadata.qs2, so re-running on 5/26
WTmerged.obj <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260421_mergedobj_no_metadata.qs2") 

# 02. attach metadata ----
## 02a. read in metadata excel sheet ----
WT_metadata <- data.frame(
  read_excel("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/WT_AgingXMetabolism_Xenium_metadata.xlsx")) 

## 02b. extract metadata from merged seurat obj ----
meta <- WTmerged.obj@meta.data

## 02c. convert rownames to cell id ----
meta <- meta %>%
  rownames_to_column(var = "cell_id")

## 02d. join metadata excel sheet to merged obj metadata ----
meta <- meta %>%
  left_join(WT_metadata, by = "sample_ID"
            #, keep = T
  )

## 02e. revert cell id back to columns ----
meta <- column_to_rownames(meta, var = "cell_id")

## 02f. attach metadata back to merged obj ----
WTmerged.obj@meta.data <- meta


WTmerged.obj <- UpdateSeuratObject(WTmerged.obj)


# 03. save merged obj not processed but w/ metadata attached ----
qs_save(WTmerged.obj, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260526_mergedobj_w_metadata_noprocessing.qs2")


# 04. split layers by slide ----
WTmerged.obj[["Xenium"]] <- split(WTmerged.obj[["Xenium"]], f = WTmerged.obj$slide)

print(WTmerged.obj)

# 05. run SCT and PCA ----
WTmerged.obj <- SCTransform(WTmerged.obj, assay = "Xenium")

WTmerged.obj <- RunPCA(WTmerged.obj, npcs = 40)


# 06. Integrate layers using harmony method  ----
# following seurat vignette https://satijalab.org/seurat/articles/seurat5_integration 

WTmerged.obj <- IntegrateLayers(object = WTmerged.obj, 
                                method = HarmonyIntegration, 
                                orig.reduction = "pca", 
                                new.reduction = "harmony")

# 07. save harmony obj without the rest of the processing done on it ----
qs_save(WTmerged.obj, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260526_harmonyintlayers_noprocessing.qs2")

# 08. standard seurat pipeline ----

WTmerged.obj <- FindNeighbors(WTmerged.obj, reduction = "harmony", dims = 1:30)

WTmerged.obj <- FindClusters(WTmerged.obj, resolution = 0.5) 

WTmerged.obj <- RunUMAP(WTmerged.obj, reduction = "harmony", dims = 1:30)

# 05. save fully processed harmony integrated obj ----
qs_save(WTmerged.obj, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260526_harmonyint_fullyprocessed.qs2")
