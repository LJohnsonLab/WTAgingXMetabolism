########################################################################
# Name: CCL_ST_SpaNorm_RPCA_WTAgingXMet_13.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 13. astrocyte subclusters metabolic pathway figures 
# Input Files: 
# Final Output Files: 
# Date created: 6/23/26
# Last updated: 7/7/26
# Author: Chloe Lucido
########################################################################


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
library(pheatmap) # for pretty heatmaps 

# set seed
set.seed(42)


######### Paths ############
zscore_metadata_inputobj <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/SpaNorm_RPCA/metabolic_scored_metadata/20260623_SpaNorm_RPCA_astro_zscore_metadata.qs2"
fig_path <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/heatmaps/"
date <- "20260623_" 
###########################

# 01. read in metadata objs ----
zscore_metadata <- qs_read(zscore_metadata_inputobj)

## 01a. extract cols ----
zscore_cols <- grep("_zscore", colnames(zscore_metadata), value = TRUE)
global_zscore_cols <- grep("_global_zscore", zscore_cols, value = TRUE)
celltype_zscore_cols <- grep("_celltype_zscore", zscore_cols, value = TRUE)


## 01b. make age group color palette ----
# age group palette 
age_cols_palette <- c(
  "16" = "#F2AF4A", 
  "36" = "#EB7F54", 
  "59" = "#C36377", 
  "92" = "#61599D"
)

# 02. GLOBAL Z-SCORE HEATMAPs ----
########################################################################################
## by cell type only ----


# global heatmap matrix
global_zscore_cell_heatmap_mat <- zscore_metadata %>%
  select(astro_subtype, all_of(global_zscore_cols)) %>%
  group_by(astro_subtype) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames("astro_subtype") %>%
  t() %>%  # transpose
  as.matrix()

# 3. Clean up row names
rownames(global_zscore_cell_heatmap_mat) <- gsub("_global_zscore", "", rownames(global_zscore_cell_heatmap_mat))

# 4. Symmetric color breaks so 0 = white
global_zscore_cell_max_abs <- max(abs(global_zscore_cell_heatmap_mat), na.rm = TRUE)

global_zscore_bycell_heatmap <- pheatmap(
  global_zscore_cell_heatmap_mat,
  color        = colorRampPalette(c("blue", "white", "red"))(100),
  breaks       = seq(-global_zscore_cell_max_abs, global_zscore_cell_max_abs, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  angle_col = 45,
  fontsize_row = 10,
  fontsize_col = 10,
  main         = "Metabolic pathway global z-scores by cell type"
)



########################################################################################
## by age only ----
# global heatmap matrix
global_zscore_age_heatmap_mat <- zscore_metadata %>%
  select(Age, all_of(global_zscore_cols)) %>%
  group_by(Age) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames("Age") %>%
  t() %>%  # tranpose
  as.matrix()


# 3. Clean up row names
rownames(global_zscore_age_heatmap_mat) <- gsub("_global_zscore", "", rownames(global_zscore_age_heatmap_mat))

# 4. Symmetric color breaks so 0 = white
global_zscore_age_max_abs <- max(abs(global_zscore_age_heatmap_mat), na.rm = TRUE)

global_zscore_byage_heatmap <- pheatmap(
  global_zscore_age_heatmap_mat,
  color        = colorRampPalette(c("blue", "white", "red"))(100),
  breaks       = seq(-global_zscore_age_max_abs, global_zscore_age_max_abs, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,  # keep age groups in order
  border_color = NA,
  angle_col = 315,
  fontsize_row = 10,
  fontsize_col = 10,
  main         = "Metabolic pathway global z-scores by age group"
)


########################################################################################
## by age and cell type ----

# global heatmap matrix
global_zscore_agecell_heatmap_mat <- zscore_metadata %>%
  select(astro_subtype, Age, all_of(global_zscore_cols)) %>%
  group_by(astro_subtype, Age) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  unite("group", astro_subtype, Age, sep = "_") %>%  # combine into one label
  column_to_rownames("group") %>%
  t() %>%
  as.matrix()

# 3. Clean up row names
rownames(global_zscore_agecell_heatmap_mat) <- gsub("_global_zscore", "", rownames(global_zscore_agecell_heatmap_mat))

# 3. Annotation bar — split on known age suffixes FOR AGE AND CELL TYPE
zscore_annotation_col <- data.frame(
  #Cell_Type = gsub("_(16|36|59|92)$", "", colnames(zscore_agecell_heatmap_mat)),
  Age_Group = gsub(".*_(16|36|59|92)$", "\\1", colnames(global_zscore_agecell_heatmap_mat)),
  row.names = colnames(global_zscore_agecell_heatmap_mat)
)

# 4. Symmetric color breaks so 0 = white
global_zscore_agecell_max_abs <- max(abs(global_zscore_agecell_heatmap_mat), na.rm = TRUE)

global_zscore_byagecell_heatmap <- pheatmap(
  global_zscore_agecell_heatmap_mat,
  color        = colorRampPalette(c("blue", "white", "red"))(100),
  breaks       = seq(-global_zscore_agecell_max_abs, global_zscore_agecell_max_abs, length.out = 101),
  annotation_col = zscore_annotation_col,
  annotation_colors = list(
    Age_Group = age_cols_palette   # name must match column name in annotation_col exactly
  ),
  angle_col = 45,
  cluster_rows = FALSE,
  cluster_cols = FALSE,  # keep age groups in order
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10,
  main         = "Metabolic pathway global z-scores by cell type and age"
)


# 03. CELL INDEPENDENT Z-SCORE HEATMAPs ----
########################################################################################
## cell independent zscore by cell type only ----

# cell-specific heatmap matrix
cellind_zscore_cell_heatmap_mat <- zscore_metadata %>%
  select(astro_subtype, all_of(celltype_zscore_cols)) %>%
  group_by(astro_subtype) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames("astro_subtype") %>%
  t() %>%  # transpose
  as.matrix()

# 3. Clean up row names
rownames(cellind_zscore_cell_heatmap_mat) <- gsub("_celltype_zscore", "", rownames(cellind_zscore_cell_heatmap_mat))

# 4. Symmetric color breaks so 0 = white
cellind_zscore_cell_max_abs <- max(abs(cellind_zscore_cell_heatmap_mat), na.rm = TRUE)

cellind_zscore_bycell_heatmap <- pheatmap(
  cellind_zscore_cell_heatmap_mat,
  color        = colorRampPalette(c("blue", "white", "red"))(100),
  breaks       = seq(-cellind_zscore_cell_max_abs, cellind_zscore_cell_max_abs, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  angle_col = 45,
  fontsize_row = 10,
  fontsize_col = 10,
  main         = "Metabolic pathway z-scores by cell type"
)


########################################################################################
## cell independent zscore by age only ----
# global heatmap matrix
cellind_zscore_age_heatmap_mat <- zscore_metadata %>%
  select(Age, all_of(celltype_zscore_cols)) %>%
  group_by(Age) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames("Age") %>%
  t() %>%  # tranpose
  as.matrix()


# 3. Clean up row names
rownames(cellind_zscore_age_heatmap_mat) <- gsub("_celltype_zscore", "", rownames(cellind_zscore_age_heatmap_mat))

# 4. Symmetric color breaks so 0 = white
cellind_zscore_age_max_abs <- max(abs(cellind_zscore_age_heatmap_mat), na.rm = TRUE)

cellind_zscore_byage_heatmap <- pheatmap(
  cellind_zscore_age_heatmap_mat,
  color        = colorRampPalette(c("blue", "white", "red"))(100),
  breaks       = seq(-cellind_zscore_age_max_abs, cellind_zscore_age_max_abs, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,  # keep age groups in order
  border_color = NA,
  angle_col = 45,
  fontsize_row = 10,
  fontsize_col = 10,
  main         = "Metabolic pathway z-scores by age group"
)



########################################################################################
## cell independent zscore by age and cell type ----

# cell-indepedent heatmap matrix
cellind_zscore_agecell_heatmap_mat <- zscore_metadata %>%
  select(astro_subtype, Age, all_of(celltype_zscore_cols)) %>%
  group_by(astro_subtype, Age) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  unite("group", astro_subtype, Age, sep = "_") %>%  # combine into one label
  column_to_rownames("group") %>%
  t() %>%
  as.matrix()

# 3. Clean up row names
rownames(cellind_zscore_agecell_heatmap_mat) <- gsub("_celltype_zscore", "", rownames(cellind_zscore_agecell_heatmap_mat))

# 3. Annotation bar — split on known age suffixes FOR AGE AND CELL TYPE
zscore_annotation_col <- data.frame(
  #Cell_Type = gsub("_(16|36|59|92)$", "", colnames(zscore_agecell_heatmap_mat)),
  Age_Group = gsub(".*_(16|36|59|92)$", "\\1", colnames(global_zscore_agecell_heatmap_mat)),
  row.names = colnames(cellind_zscore_agecell_heatmap_mat)
)

# 4. Symmetric color breaks so 0 = white
cellind_zscore_agecell_max_abs <- max(abs(cellind_zscore_agecell_heatmap_mat), na.rm = TRUE)

cellind_zscore_byagecell_heatmap <- pheatmap(
  cellind_zscore_agecell_heatmap_mat,
  color        = colorRampPalette(c("blue", "white", "red"))(100),
  breaks       = seq(-cellind_zscore_agecell_max_abs, cellind_zscore_agecell_max_abs, length.out = 101),
  annotation_col = zscore_annotation_col,
  annotation_colors = list(
    Age_Group = age_cols_palette   # name must match column name in annotation_col exactly
  ),
  angle_col = 45,
  cluster_rows = FALSE,
  cluster_cols = FALSE,  # keep age groups in order
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10,
  main         = "Metabolic pathway z-scores by cell type and age"
)


# save heatmaps ----
## global z score heatmaps ----
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_global_zscore_bycell_heatmap.png"), plot = global_zscore_bycell_heatmap, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_global_zscore_byage_heatmap.png"), plot = global_zscore_byage_heatmap, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_global_zscore_byagecell_heatmap.png"), plot = global_zscore_byagecell_heatmap, path = fig_path, height = 8, width = 11, dpi = 600)

## cell-specific z score heatmaps ----
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_cellind_zscore_bycell_heatmap.png"), plot = cellind_zscore_bycell_heatmap, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_cellind_zscore_byage_heatmap.png"), plot = cellind_zscore_byage_heatmap, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_cellind_zscore_byagecell_heatmap.png"), plot = cellind_zscore_byagecell_heatmap, path = fig_path, height = 8, width = 11, dpi = 600)



