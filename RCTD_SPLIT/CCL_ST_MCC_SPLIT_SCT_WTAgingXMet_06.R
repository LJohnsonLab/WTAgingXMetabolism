########################################################################
# Name: CCL_ST_MCC_SPLIT_SCT_WTAgingXMet_06.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 06. clean (remove cells <5 nCount), attach images to split shift obj, and prep for refRPCA int (run SCT per slide) 
# Input:  
# Final output: 
# Date created: 8/7/26
# Last updated: 8/7/26
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

# Use cores requested on SLURM job script
plan("multicore", workers = 8)
options(future.globals.maxSize = 400 * 1024^3)  # slightly under your --mem

# to prevent this error: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so.
options(future.rng.onMisuse = "resolve")
options(future.seed = TRUE)

# set seed
set.seed(42)

###### PATHS ########
split_shift_obj <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260805_SPLIT_shift.qs2"
source_obj <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260515_mergedobj_w_metadata_noprocessing.qs2"
output_dir <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/RCTD_SPLIT_refRPCA/"
date <- "20260807_"
#####################

split.obj <- qs_read(split_shift_obj)

source.obj <- qs_read(source_obj)

# 00. clean obj ----
split.obj <- subset(split.obj, subset = nCount_Xenium > 5)

dim(split.obj)
# 474 1522573

# 01. attach images back to split shift obj ----

## 01a. extract fov names from source obj (source.obj)
fov_names <- Images(source.obj)

## 01b. loop through fovs in sources obj and subset to common cell ids btwn source and split obj
fovs_matched <- imap(set_names(fov_names), function(fov, fov_name) {
  
  ### identify common cells between source obj (source.obj) and split obj
  cells_in_fov <- intersect(Cells(source.obj[[fov_name]]), colnames(split.obj))
  
  ### returning NULL for fovs that have no common cell ids
  if (length(cells_in_fov) == 0) return(NULL)
  
  ### subsetting fovs to only cells in split.obj
  subset(source.obj[[fov_name]], cells = cells_in_fov)
}) |> 
  compact() # drop any FOVs with no cells

## 01c. SAVE fovs_matched 
qs_save(fovs_matched, paste0(output_dir, date, "fovs_matched_06.qs2"))

## 01d. attach matched fovs to split.obj
for(fov_name in names(fovs_matched)) {
  split.obj[[fov_name]] <- fovs_matched[[fov_name]]
}

### SANITY CHECK: confirm images attached 
Images(split.obj)

## 01e. SAVE split obj with images attached
qs_save(split.obj, paste0(output_dir, date, "split_shift_w_imgs.qs2"))

  
# 02. prepare for ref RPCA integration ----
## 02a. Update merged obj 
split.obj <- UpdateSeuratObject(split.obj)

dim(split.obj)

## 02b. split split.obj by slide ----
obj_list <- SplitObject(split.obj, split.by = "slide_ID")

## 02c. SCTransform each slide independently ----
obj_list <- map(obj_list, ~ SCTransform(.x, assay = "Xenium",  vst.flavor = "v2", conserve.memory = TRUE, verbose = FALSE))

# 03. Save split SCTransformed obj_list ----
qs_save(obj_list, paste0(output_dir, date, "split_SCT_objlist_06.qs2"))
