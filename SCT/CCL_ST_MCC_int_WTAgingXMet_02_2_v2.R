########################################################################
# Name: CCL_ST_MCC_int_WTAgingXMet_02_2_v2.R **ONLY TO BE RUN ON SUPERCOMPUTER (MCC)**
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - MCC 02. integrating and standard seurat processing using MCC 
#             - 01 read in SCTransformed obj list 
#             - 02 prepare for RPCA integration
#             - 03 integrate
#             - 04 save integrated obj (no dim reductions, only SCT) => 20260508_nodimreduc_integratedobj_02_2.qs2
#             - 05 dim reductions (RunPCA, RunUMAP, FindNeighbors and clusters)
#             - 06 save final processed integrated obj => 20260508_final_integratedobj_02_2.qs2
# Input Files: int part 2_1 obj list with SCT =>  20260508_objlist_SCT_02_1.qs2
# Output Files: integrated processed obj => 20260508_final_integratedobj_02_2.qs2
# Date created: 5/15/26
# Last updated: 5/26/26
# Author: Chloe Lucido
# Version history:
#       - v2 uses 2 reference slides to reduce matrix size to avoid hitting 32bit-integer seurat limit, also only using 2000 features instead of 3000
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
plan("multicore", workers = 8)
options(future.globals.maxSize = 850 * 1024^3)  # slightly under your --mem

# set seed
set.seed(42)

# to prevent this error: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so.
options(future.rng.onMisuse = "resolve")
options(future.seed = TRUE)


# 01. read in objlist_SCT obj from 02_1 script ----
#obj_list <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260508_objlist_SCT_02_1.qs2")


# 02. Prepare for Integration ----
## 02a. Select integration features ----
#features <- SelectIntegrationFeatures(
#  object.list = obj_list,
#  nfeatures = 2000 
#)

## 02b. Prep SCT for integration ----
#obj_list <- PrepSCTIntegration(
#  object.list = obj_list,
#  anchor.features = features
#)

## 02c. Run PCA on each object ----
#obj_list <- map(obj_list, ~ RunPCA(.x, features = features, verbose = FALSE))

## 02d. Find integration anchors using RPCA
#anchors <- FindIntegrationAnchors(
#  object.list = obj_list,
#  normalization.method = "SCT",
#  anchor.features = features,
#  reduction = "rpca",
#  reference = c(2, 7), # using slides 0022573 and 0069118 because they contain the most cells (more anchors)
#  dims = 1:30,
#  k.anchor = 5
#)


## 2e. Save integration anchors ----
#qs_save(anchors, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260515_refRPCAint_anchors_02_2_v2.qs2")

# 03. Integrate ----
#integrated.obj <- IntegrateData(
#  anchorset = anchors,
#  normalization.method = "SCT",
#  dims = 1:30
#)

# 04. save integrated obj (no dim reductions) ----
#qs_save(integrated.obj, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260515_nodimreduc_refRPCA_integratedobj_02_2_v2.qs2") 

# pick up where I left off 
integrated.obj <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260515_nodimreduc_refRPCA_integratedobj_02_2_v2.qs2")

# 05. Processing ----

## 05a. set default assay to integrated ----

DefaultAssay(integrated.obj) <- "integrated"

## 05b. Standard Seurat pipeline ----

integrated.obj <- RunPCA(integrated.obj, npcs = 50, features = rownames(integrated.obj))

integrated.obj <- RunUMAP(integrated.obj, dims = 1:30)

integrated.obj <- FindNeighbors(integrated.obj, reduction = "pca", dims = 1:30)

### print out available graphs 
print(integrated.obj@graphs)

integrated.obj <- FindClusters(integrated.obj, resolution = 0.5)

# 06. Save SCTransformed integrated obj ----
qs_save(integrated.obj, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260526_final_refRPCA_integratedobj_02_2_v2.qs2") 
