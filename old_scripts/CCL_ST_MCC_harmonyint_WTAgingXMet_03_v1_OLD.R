########################################################################
# Name: CCL_ST_MCC_harmonyint_WTAgingXMet_03.R **ONLY TO BE RUN ON SUPERCOMPUTER (MCC)**
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
# Last updated: 5/15/26
# Author: Chloe Lucido
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
options(future.globals.maxSize = 400 * 1024^3)  # slightly under your --mem


# set seed
set.seed(42)

# to prevent this error: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so.
options(future.rng.onMisuse = "resolve")
options(future.seed = TRUE)
 

 
# 01. read in merged obj (not processed) w/ metadata attached ----
 
WTmerged.obj <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260421_mergedobj_SCT.qs2") # this obj already has pca reduction

# 02. Run Harmony integration  ----
# NOTE: 
# tried running rpca integration but the sparse matrix indices exceeded the 32-bit integer limit (~2.1 billion elements)
# rpca also used 700GB of memory and took over a day to run before failing because of the above comment 
# deciding to use Harmony because it is faster, works better for many samples/batches, and memory efficient

integrated.obj <- RunHarmony(
                             WTmerged.obj, 
                             group.by.vars = "slide", 
                             reduction = "pca", 
                             reduction.save = "harmony"
                             )

# 03. save harmony obj without the rest of the processing done on it ----
qs_save(integrated.obj, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260513_harmonyint_noprocessing.qs2")

integrated.obj <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260513_harmonyint_noprocessing.qs2")

# 04. standard seurat pipeline ----
integrated.obj <- RunUMAP(integrated.obj, reduction = "harmony", dims = 1:30)

integrated.obj <- FindNeighbors(integrated.obj, reduction = "harmony", dims = 1:30)

print(names(integrated.obj@graphs))

integrated.obj <- FindClusters(integrated.obj, resolution = 0.5, graph.name = "harmony_snn") 

# 05. save fully processed harmony integrated obj ----
qs_save(integrated.obj, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260513_harmonyint.qs2")
