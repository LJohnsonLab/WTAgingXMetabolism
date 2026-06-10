########################################################################
# Name: CCL_ST_MCC_SpaNormprocessing_noint_WTAgingXMet_03.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 03. processing SpaNorm obj (FindVariableFeatures->ScaleData) NO INTEGRATION 
# Input Files: merged obj with SpaNorm assay from   
# Final Output Files: 
# Date created: 6/04/26
# Last updated: 6/04/26
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
options(future.globals.maxSize = 850 * 1024^3)  # slightly under your --mem
 
# set seed
set.seed(42)
 
# to prevent this error: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so.
options(future.rng.onMisuse = "resolve")
options(future.seed = TRUE)


# read in SpaNormalized obj
WTmerged.obj <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/20260604_mergedobj_SpaNorm_Assay_noscaling.qs2")

DefaultAssay(WTmerged.obj) <- "SpaNorm"

# Use all 480 genes as variable features; on a targeted Xenium panel every gene was chosen as informative by design,
# and FindVariableFeatures("vst") requires a counts layer the SpaNorm assay does not carry.
VariableFeatures(WTmerged.obj) <- rownames(WTmerged.obj)

# RunPCA expects a scale.data layer; SCT writes it implicitly, but our custom SpaNorm assay only carries the data layer.
WTmerged.obj <- ScaleData(WTmerged.obj, verbose = FALSE)

# Default npcs (50) is fine; we only use the first 30 downstream.
WTmerged.obj <- RunPCA(WTmerged.obj, npcs = 50, verbose = FALSE)

# 1:30 PCs matches the SCT twins; FindNeighbors default k = 20.
WTmerged.obj <- FindNeighbors(WTmerged.obj, reduction = "pca", dims = 1:30, verbose = FALSE)

# algorithm = 1 is plain Louvain; resolution 0.5 matches the SCT twins for a like-for-like comparison.
WTmerged.obj <- FindClusters(WTmerged.obj, resolution = 0.5,  verbose = FALSE)

# UMAP only for visual QC of the clusters; not used by any downstream stat.
WTmerged.obj <- RunUMAP(WTmerged.obj, reduction = "pca", dims = 1:30, verbose = FALSE)


# attach these clusters to metadata ----
Idents(WTmerged.obj) <- "seurat_clusters"

WTmerged.obj <- AddMetaData(WTmerged.obj, metadata = Idents(WTmerged.obj), col.name = "clusters_res05")

# FindAllMarkers
SpaNorm.markers <- FindAllMarkers(WTmerged.obj, group.by = "seurat_clusters", only.pos = T, logfc.threshold = 0.25)

## save all markers as csv file ----
write_csv(SpaNorm.markers, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260604_SpaNorm_noint_res05_allmarkers.csv")

## 07b. find top 10 markers -----
SpaNorm.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  ungroup() -> top10_SpaNorm

## save top10 markers as csv ----
write_csv(top10_SpaNorm, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260604_SpaNorm_noint_res05_top10_markers.csv")

# re-save obj with cluster information ----
qs_save(WTmerged.obj, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260604_SpaNorm_noint_processed_wclus.qs2")


