########################################################################
# Name: CCL_ST_MCC_int_WTAgingXMet_02_1.R **ONLY TO BE RUN ON SUPERCOMPUTER (MCC)**
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - MCC 02. integrating and running SCTransform using MCC 
#             - 01 read in merged obj (w/ metadata, not processed)
#             - 02 update merged obj
#             - 03 prepare for RPCA integration
#             - 04 integrate data 
#             - 05 save integrated obj => 20260506_integrated_mergedobj.qs2
#             - 06 process obj (SCT)
#             - 07 save integrated and processed obj => 20260506_integratedobj_SCT.qs2
# Input Files: 20260420_mergedobjv1_01.qs2 (merged, w/ metadata, not processed)
# Output Files: integrated SCTransformed merged obj => 20260506_integratedobj_SCT.qs2
# Date created: 5/6/26
# Last updated: 5/8/26
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
plan("multicore", workers = 8)
options(future.globals.maxSize = 400 * 1024^3)  # slightly under your --mem

# set seed
set.seed(42)

# to prevent this error: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so.
options(future.rng.onMisuse = "resolve")
options(future.seed = TRUE)



# 01. read in merged obj (not processed) w/ metadata attached ----

WTmerged.obj <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/20260420_mergedobjv1_01.qs2")


# 02. Update merged obj ----

WTmerged.obj <- UpdateSeuratObject(WTmerged.obj)


# 03. Prepare obj for RPCA Integration ----

## 03a. split merged obj by slide ----
obj_list <- SplitObject(WTmerged.obj, split.by = "slide_ID")

## 03b. SCTransform each slide independently ----
obj_list <- map(obj_list, ~ SCTransform(.x, assay = "Xenium",  vst.flavor = "v2", conserve.memory = TRUE, verbose = FALSE))


# 04. Save SCTransformed obj_list ----
qs_save(obj_list, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260508_objlist_SCT_02_1.qs2")
