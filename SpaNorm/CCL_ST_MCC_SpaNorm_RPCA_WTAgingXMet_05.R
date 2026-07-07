########################################################################
# Name: CCL_ST_SpaNorm_RPCA_WTAgingXMet_05.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - SpaNorm 05. subclustering clusters 9, 18, 1 and 3 
#             - 01 
#             - 02 
#             - 03 
#             - 04 
#             - 05 
#             - 06  
#             - 07 
#             - 08  
# Input Files:   
# Final Output Files: 
# Date created: 6/9/26
# Last updated: 7/7/26
# Author: Chloe Lucido
########################################################################

# Auto-install missing packages
required_packages <- c("Seurat", "future", "ggplot2", "readr", "dplyr", 
                       "data.table","tidyr", "tibble", "stringr", "forcats", "lubridate",  "readxl", "patchwork", 
                       "qs2", "RColorBrewer", "Polychrome", "purrr",
                       "presto", "glmGamPoi", "ragg")

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
library(ragg) # for saving images on mcc

# ONLY TO BE RUN ON MCC 

# Use the cores you requested
plan("multicore", workers = 8)
options(future.globals.maxSize = 400 * 1024^3)  # slightly under your --mem

# set seed
set.seed(42)

# to prevent this error: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so.
options(future.rng.onMisuse = "resolve")
options(future.seed = TRUE)

# set ggsave to use ragg globally in order to save figures 
options(ragg.use_agg = TRUE)


#########Paths##########
root_path_to_obj <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260605_SpaNorm_RPCAint_processed.qs2"
output_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/"
fig_output_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/SpaNorm_RPCA_figs/"
csv_output_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/SpaNorm_cluster_markers/"
date <- "20260610_"
#######################


# read in obj ----
SpaNorm_RPCA.obj <- qs_read(root_path_to_obj)

Idents(SpaNorm_RPCA.obj) <- "clusters_res05"

# print out graph names 
print(SpaNorm_RPCA.obj@graphs)

#output:
#$integrated_nn
#A Graph object containing 1609408 cells
#$integrated_snn
#A Graph object containing 1609408 cells

# Find sub clusters for cluster 9
# NOTE: cluster 9 has VLMC, BAM and fibroblast markers in top10, also hits for astro markers (Gfap, Aldh1l1, and Aqp4)
SpaNorm_RPCA.obj <- FindSubCluster(SpaNorm_RPCA.obj, cluster = "9", resolution = 0.2, graph.name = "integrated_snn")

print(levels(Idents(SpaNorm_RPCA.obj)))

print(colnames(SpaNorm_RPCA.obj@meta.data))

Idents(SpaNorm_RPCA.obj) <- "sub.cluster"

print(levels(Idents(SpaNorm_RPCA.obj)))

SpaNorm_RPCA.obj$clusters_sub9 <- Idents(SpaNorm_RPCA.obj)

Idents(SpaNorm_RPCA.obj) <- "clusters_sub9"

# find subclusters for cluster 18
# NOTE: cluster 18 has multiple neuronal-type markers and oligo markers, but clusters with oligos on UMAP 
SpaNorm_RPCA.obj <- FindSubCluster(SpaNorm_RPCA.obj, cluster = "18", resolution = 0.2, graph.name = "integrated_snn")

Idents(SpaNorm_RPCA.obj) <- "sub.cluster"

print(levels(Idents(SpaNorm_RPCA.obj)))

SpaNorm_RPCA.obj$clusters_sub9_18 <- Idents(SpaNorm_RPCA.obj)

Idents(SpaNorm_RPCA.obj) <- "clusters_sub9_18"

# need to joinlayers for FindAllMarkers to work
SpaNorm_RPCA.obj <- JoinLayers(SpaNorm_RPCA.obj)

# find markers

SpaNorm_RPCA.submarkers <- FindAllMarkers(SpaNorm_RPCA.obj, group.by = "clusters_sub9_18")

print(colnames(SpaNorm_RPCA.submarkers))

SpaNorm_RPCA.submarkers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  ungroup() -> top10_sub

write_csv(top10_sub, paste0(csv_output_path, date, "top10_SpaNormRPCA_sub9_18.csv"))

# figures 
UMAP_sub <- DimPlot(SpaNorm_RPCA.obj, group.by = "clusters_sub9_18", cols = "polychrome")
imgdim_sub <- ImageDimPlot(SpaNorm_RPCA.obj, split.by = "clusters_sub9_18", cols = "polychrome")


# save figures
ggsave(filename = paste0(date, "SpaNormRPCA_UMAP_sub9_18.png"), plot = UMAP_sub, path = fig_output_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNormRPCA_imgdim_sub9_18.png"), plot = imgdim_sub, path = fig_output_path, width = 11, height = 10)

# save object with subcluster info for clusters 9 and 18
qs_save(SpaNorm_RPCA.obj, paste0(output_path, date, "SpaNorm_RPCAint_processed.qs2"))

# check cell counts across subclusters 
SpaNorm_RPCA.obj <- qs_read(paste0(output_path, "20260609_SpaNorm_RPCAint_processed.qs2"))

print(table(SpaNorm_RPCA.obj$clusters_sub9_18))

print(table(SpaNorm_RPCA.obj$clusters_sub9_18, SpaNorm_RPCA.obj$sample_ID))

# FURTHER SUBCLUSTERING
## NOTES: 
## subclustering 18_0 and 18_1 => only half of these clusters hit for canonical oligo markers (based on feature plots)
## subclustering 1 to try and get rid of the cells that are clustering with cluster 2 => these cells are hitting heavily for both canonical oligo and endo markers (clearly segmentation artifacts)
## subclustering 3 (astrocytes) because some cells are clustering with 18_0 and 18_1 (bad oligo clusters), those same cells aren't hitting heavily with the canonical astro markers 

# subclustering 18_0
Idents(SpaNorm_RPCA.obj) <- "clusters_sub9_18"

SpaNorm_RPCA.obj <- FindSubCluster(SpaNorm_RPCA.obj, cluster = "18_0", resolution = 0.1, graph.name = "integrated_snn")

Idents(SpaNorm_RPCA.obj) <- "sub.cluster"

SpaNorm_RPCA.obj$clusters_sub18_0 <- Idents(SpaNorm_RPCA.obj)

Idents(SpaNorm_RPCA.obj) <- "clusters_sub18_0"


# subclustering 18_1
SpaNorm_RPCA.obj <- FindSubCluster(SpaNorm_RPCA.obj, cluster = "18_1", resolution = 0.1, graph.name = "integrated_snn")

Idents(SpaNorm_RPCA.obj) <- "sub.cluster"

SpaNorm_RPCA.obj$clusters_sub18_0and1 <- Idents(SpaNorm_RPCA.obj)

Idents(SpaNorm_RPCA.obj) <- "clusters_sub18_0and1"


# subclustering 1 
SpaNorm_RPCA.obj <- FindSubCluster(SpaNorm_RPCA.obj, cluster = "1", resolution = 0.2, graph.name = "integrated_snn")

Idents(SpaNorm_RPCA.obj) <- "sub.cluster"

SpaNorm_RPCA.obj$clusters_sub1 <- Idents(SpaNorm_RPCA.obj)

Idents(SpaNorm_RPCA.obj) <- "clusters_sub1"

# subclustering 3
SpaNorm_RPCA.obj <- FindSubCluster(SpaNorm_RPCA.obj, cluster = "3", resolution = 0.2, graph.name = "integrated_snn")

Idents(SpaNorm_RPCA.obj) <- "sub.cluster"

SpaNorm_RPCA.obj$clusters_sub3 <- Idents(SpaNorm_RPCA.obj)

Idents(SpaNorm_RPCA.obj) <- "clusters_sub3"

# Find Markers 

# need to joinlayers for FindAllMarkers to work
SpaNorm_RPCA.obj <- JoinLayers(SpaNorm_RPCA.obj)

SpaNorm_RPCA.submarkers <- FindAllMarkers(SpaNorm_RPCA.obj, group.by = "clusters_sub3")

print(colnames(SpaNorm_RPCA.submarkers))

SpaNorm_RPCA.submarkers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  ungroup() -> top10_sub

write_csv(top10_sub, paste0(csv_output_path, date, "top10_SpaNormRPCA_sub18_0_18_1_1_3.csv"))

# figures 
UMAP_sub <- DimPlot(SpaNorm_RPCA.obj, group.by = "clusters_sub3", cols = "polychrome", label = T)
imgdim_split_sub <- ImageDimPlot(SpaNorm_RPCA.obj, split.by = "clusters_sub3", cols = "polychrome")
imgdim_sub <- ImageDimPlot(SpaNorm_RPCA.obj, group.by = "clusters_sub3", cols = "polychrome")

# save figures 
ggsave(filename = paste0(date, "SpaNormRPCA_UMAP_sub18_0_18_1_1_3.png"), plot = UMAP_sub, path = fig_output_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNormRPCA_imgdim_sub18_0_18_1_1_3.png"), plot = imgdim_sub, path = fig_output_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNormRPCA_imgdim_split_sub18_0_18_1_1_3.png"), plot = imgdim_split_sub, path = fig_output_path, width = 11, height = 10)

# print cell counts per cluster 
print(table(SpaNorm_RPCA.obj$clusters_sub3))

# save obj
qs_save(SpaNorm_RPCA.obj, paste0(output_path, date, "SpaNorm_RPCAint_processed_subclustered.qs2"))

# read in last obj
SpaNorm_RPCA.obj <- qs_read(paste0(output_path, date, "SpaNorm_RPCAint_processed_subclustered.qs2"))

UMAP_sub <- DimPlot(SpaNorm_RPCA.obj, group.by = "clusters_sub3")
imgdim_split_sub <- ImageDimPlot(SpaNorm_RPCA.obj, split.by = "clusters_sub3")
imgdim_sub <- ImageDimPlot(SpaNorm_RPCA.obj, group.by = "clusters_sub3")

# save figures
ggsave(filename = paste0(date, "SpaNormRPCA_UMAP_sub18_0_18_1_1_3.png"), plot = UMAP_sub, path = fig_output_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNormRPCA_imgdim_sub18_0_18_1_1_3.png"), plot = imgdim_sub, path = fig_output_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNormRPCA_imgdim_split_sub18_0_18_1_1_3.png"), plot = imgdim_split_sub, path = fig_output_path, width = 11, height = 10)


