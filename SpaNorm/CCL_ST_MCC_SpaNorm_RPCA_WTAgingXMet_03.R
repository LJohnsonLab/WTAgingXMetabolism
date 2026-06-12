########################################################################
# Name: CCL_ST_MCC_SpaNorm_RPCAint_WTAgingXMet_03.R **ONLY TO BE RUN ON SUPERCOMPUTER (MCC)**
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - MCC 02. running RPCA integration on SpaNorm obj 
# Input Files: 
# Output Files: 
# Date created: 6/5/26
# Last updated: 6/5/26
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
plan("sequential")
options(future.globals.maxSize = 850 * 1024^3)  # slightly under your --mem

# set seed
set.seed(42)

# to prevent this error: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so.
options(future.rng.onMisuse = "resolve")
options(future.seed = TRUE)


# 01. read in objlist_SCT obj from 02_1 script ----
SpaNorm_merged.obj <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/20260604_mergedobj_SpaNorm_Assay_noscaling.qs2")

# 02. Split merged obj by sample_ID into a list ----
DefaultAssay(SpaNorm_merged.obj) <- "SpaNorm"

SpaNorm_obj_list <- SplitObject(SpaNorm_merged.obj, split.by = "sample_ID")

# 03. Prepare for Integration ----
## 03a. Select integration features ----
features <- SelectIntegrationFeatures(
  object.list = SpaNorm_obj_list,
  nfeatures = 2000
)


## 03b. Run PCA on each sample ----
SpaNorm_obj_list <- map(SpaNorm_obj_list, ~ {
  .x <- FindVariableFeatures(.x)
  .x <- ScaleData(.x, features = features)
  .x <- RunPCA(.x, features = features, verbose = FALSE)
  .x
})

# check sample order to confirm reference indices
print(names(SpaNorm_obj_list))

qs_save(SpaNorm_obj_list, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260605_processed_SpaNorm_sampleobjlist.qs2" 
# start where I left off 
SpaNorm_obj_list <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260605_processed_SpaNorm_sampleobjlist.qs2")


# forgot to save features from earlier but since this is a targeted 480 gene panel, I can just use all genes as the features
features <- rownames(SpaNorm_obj_list[[1]])

## 03d. Find integration anchors using RPCA
anchors <- FindIntegrationAnchors(
  object.list = SpaNorm_obj_list,
  normalization.method = "LogNormalize",
  anchor.features = features,
  reduction = "rpca",
  reference = c(2, 21), # using 59wk_F2 from slide 0022573 and 59wk_M5 from slide 0069118 because they contain the most cells (more anchors)
  dims = 1:30,
  k.anchor = 5
)


## 2e. Save integration anchors ----
qs_save(anchors, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260605_SpaNorm_refRPCAint_anchors.qs2")

anchors <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260605_SpaNorm_refRPCAint_anchors.qs2")

# 03. Integrate ----
SpaNorm_int.obj <- IntegrateData(
  anchorset = anchors,
  normalization.method = "LogNormalize",
  dims = 1:30
)

# 04. save integrated obj (no dim reductions) ----
qs_save(SpaNorm_int.obj, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260605_SpaNorm_RPCA_nodimreduc.qs2")

# pick up where I left off 
SpaNorm_int.obj <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260605_SpaNorm_RPCA_nodimreduc.qs2")

# 05. Processing ----

## 05a. set default assay to integrated ----

DefaultAssay(SpaNorm_int.obj) <- "integrated"

# scale again
SpaNorm_int.obj <- ScaleData(SpaNorm_int.obj, features = rownames(SpaNorm_int.obj))

## 05b. Standard Seurat pipeline ----

SpaNorm_int.obj <- RunPCA(SpaNorm_int.obj, npcs = 50, features = rownames(SpaNorm_int.obj))

SpaNorm_int.obj <- RunUMAP(SpaNorm_int.obj, dims = 1:30)

SpaNorm_int.obj <- FindNeighbors(SpaNorm_int.obj, reduction = "pca", dims = 1:30)

### print out available graphs 
print(SpaNorm_int.obj@graphs)

SpaNorm_int.obj <- FindClusters(SpaNorm_int.obj, resolution = 0.5)


# switch to SpaNorm assay for marker finding — integrated assay lacks counts layer
DefaultAssay(SpaNorm_int.obj) <- "SpaNorm"

# attach clusters to metadata 
Idents(SpaNorm_int.obj) <- "seurat_clusters"

SpaNorm_int.obj <- AddMetaData(SpaNorm_int.obj, metadata = Idents(SpaNorm_int.obj), col.name = "clusters_res05")

# 06. Save SpaNorm RPCA integrated obj ----
qs_save(SpaNorm_int.obj, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260605_SpaNorm_RPCAint_processed.qs2")

# load in previous obj ----
SpaNorm_int.obj <- qs_read("/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260605_SpaNorm_RPCAint_processed.qs2")

SpaNorm_int.obj <- JoinLayers(SpaNorm_int.obj)

# 07. FindAllMarkers 
SpaNorm_RPCA.markers <- FindAllMarkers(SpaNorm_int.obj, group.by = "clusters_res05", only.pos = T, logfc.threshold = 0.25)

## 07a. save all markers as csv file ----
write_csv(SpaNorm_RPCA.markers, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260605_SpaNorm_RPCA_res05_allmarkers.csv")

## 07b. find top 10 markers -----
SpaNorm_RPCA.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  ungroup() -> top10_SpaNorm_refRPCA

## 07c. save top10 markers as csv ----
write_csv(top10_SpaNorm_refRPCA, "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260605_SpaNorm_RPCA_res05_top10_markers.csv")

