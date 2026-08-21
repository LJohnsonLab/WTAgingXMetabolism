########################################################################
# Name: CCL_ST_SpaNorm_RPCA_WTAgingXMet_08.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 08. stats on cell proportions and cell props line graphs and met score heatmaps 
# Input Files: 
# Final Output Files: 
# Date created: 8/17/26
# Last updated: 8/17/26
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
library(speckle) # for propeller
library(pheatmap) # for pretty heatmaps 


# ONLY TO BE RUN ON MCC 
options(future.globals.maxSize = 100 * 1024^3)  # slightly under your --mem

# set seed
set.seed(42)


######### Paths ############
props_obj <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260817_immunevasc_props.qs2" 
zscore_metadata_inputobj <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260817_SpaNorm_RPCA_immunevasc_zscore_metadata.qs2"
metadata_xlsx <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/metadata/WT_AgingXMetabolism_Xenium_metadata.xlsx"
output_path <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/SpaNorm_figs/"
date <- "20260817_"
###########################

# col palette 
## col palette ----
col_palette <- c(
  "Endothelial cells" = "#FDBF6F",
  "Pericytes" = "#C45E00",
  "Astrocytes" = "#E64B35",
  "VLMC" = "#FF7F00",
  "Microglia" = "#00897B",
  "VSMC" = "#FB9A99",
  "BAM" = "#5C8A3C",
  "T cells" = "#1F78B4"
)


# CELL PROPORTIONS WITH PROPELLER----

props <- qs_read(props_obj)
metadata <- read_xlsx(metadata_xlsx)

# The raw (untransformed) proportions are stored here:
raw_props <- props$Proportions  # rows = cell types, columns = samples

# ---- Step 2: Add Age group info per sample ----
# Build a sample -> age lookup from your Seurat metadata
sample_age <- metadata %>%
  select(sample_ID, Age) %>%
  distinct()

# Reshape proportions to long format
props_long <- as.data.frame(raw_props) %>%
  rename(CellType = clusters, sample_ID = sample, Proportion = Freq) %>%
  filter(CellType != "T cells")

# create age column using sample IDs
props_long <- props_long %>%
  mutate(Age = as.numeric(gsub("wk.*", "", sample_ID)))  # extracts 16, 36, 59, 92


# ---- Step 3: Calculate Mean, SD, SE per cell type per age ----
fine_summary_df <- props_long %>%
  group_by(CellType, Age) %>%
  summarise(
    Mean = mean(Proportion),
    SD   = sd(Proportion),
    N    = n(),
    SE   = SD / sqrt(N),
    .groups = "drop"
  )


print(fine_summary_df)


# ---- Significance labels from your propeller output ----
# Add all your cell types and their FDR values here (from: /Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260817_immunevasc_props_stats.csv)
fine_sig_df <- data.frame(
  CellType = c("Microglia", "Astrocytes"),        # MANUALLY ADD
  FDR      = c(0.000146, 0.049459746)    # MANUALLY ADD
) %>%
  mutate(
    sig_label = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ "ns"
    )
  )

# ---- Join sig labels to summary_df ----
fine_summary_df <- fine_summary_df %>%
  left_join(fine_sig_df, by = "CellType")


# normalize to % of 16 weeks for better visualization ----
fine_summary_df_normalized <- fine_summary_df %>%
  group_by(CellType) %>%
  mutate(
    baseline  = Mean[Age == 16],
    Mean_norm = (Mean / baseline) * 100,
    SE_norm   = (SE / baseline) * 100
  ) %>%
  ungroup()

# plotting by cell category -----

library(extrafont)  # for Arial font
# font_import()     # run this once the first time to import system fonts
loadfonts(device = "pdf")  # or "win" if on Windows


plot_category_norm <- function(data) {
  
  line_end_labels <- data %>%
    filter(Age == 92, !is.na(sig_label), sig_label != "ns")
  
  ggplot(data, aes(x = Age, y = Mean_norm, color = CellType, fill = CellType, group = CellType)) +
    geom_hline(yintercept = 100, linetype = "dashed",
               color = "gray50", linewidth = 0.5) +
    geom_ribbon(aes(ymin = Mean_norm - SE_norm, ymax = Mean_norm + SE_norm),
                alpha = 0.15, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 1, show.legend = FALSE) +
    geom_point(size = 3, shape = 21,
               color = "white", stroke = 0.3,
               show.legend = TRUE) +
    geom_text(data = line_end_labels,
              aes(x = 95, y = Mean_norm, label = sig_label, color = CellType),
              inherit.aes = FALSE, size = 7, show.legend = FALSE) +
    scale_color_manual(values = col_palette) +
    scale_fill_manual(values = col_palette) +
    scale_x_continuous(breaks = c(16, 36, 59, 92), limits = c(16, 100)) +
    scale_y_continuous(
      breaks = function(x) unique(c(scales::breaks_pretty()(x), 100)),
      expand = expansion(mult = c(0.05, 0.1))
    ) +
    guides(
      color = guide_legend(
        title = "Cell Type",
        override.aes = list(shape = 21, 
                            size = 4,
                            fill = col_palette[sort(unique(data$CellType))],
                            color = "white", 
                            stroke = 0.3,
                            linetype = 0)
      ),
      fill = "none"    # suppress the fill legend
    ) +
    labs(
      title = "Immune-vascular cells",
      x     = "Age (weeks)",
      y     = "% of 16 weeks",
      color = "Cell Type"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      text              = element_text(family = "Arial", color = "black"),
      axis.text         = element_text(family = "Arial", color = "black", size = 11),
      axis.title        = element_text(family = "Arial", color = "black", size = 13),
      plot.title        = element_text(family = "Arial", color = "black", size = 14,
                                       face = "bold", hjust = 0.5),
      legend.text       = element_text(family = "Arial", color = "black", size = 10),
      legend.title      = element_text(family = "Arial", color = "black", size = 11),
      axis.line         = element_line(color = "black", linewidth = 0.6),
      axis.ticks        = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(3, "pt"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_blank(),
      panel.border      = element_blank(),
      legend.position   = "right",
      legend.key        = element_blank()
    )
}



# ---- Step 4: Generate each plot ----
p_immunevasc <- plot_category_norm(fine_summary_df_normalized)


# View individually
p_immunevasc # saving with geom_ribbon alpha = 0.1

ggsave(paste0(output_path, date, "norm_to_16wks_immunevasc_props_linegraph.png"), plot = p_immunevasc, width = 11, height = 8, dpi = 600)


# METABOLIC SCORING (only doing global and cell-dependent z scores)--------------------------------------------------------------------

# 01. read in metadata objs ----
zscore_metadata <- qs_read(zscore_metadata_inputobj) %>%
  filter(immune_vasc_ann != "T cells")

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


# 03. CELL INDEPENDENT Z-SCORE HEATMAPs ----
########################################################################################
## cell independent zscore by cell type only ----

# cell-specific heatmap matrix
cellind_zscore_cell_heatmap_mat <- zscore_metadata %>%
  filter(Age == "16") %>%
  select(immune_vasc_ann, all_of(celltype_zscore_cols)) %>%
  group_by(immune_vasc_ann) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames("immune_vasc_ann") %>%
  t() %>%  # transpose
  as.matrix()

# 3. Clean up row names
rownames(cellind_zscore_cell_heatmap_mat) <- gsub("_celltype_zscore", "", rownames(cellind_zscore_cell_heatmap_mat))

# 4. Symmetric color breaks so 0 = white
cellind_zscore_cell_max_abs <- max(abs(cellind_zscore_cell_heatmap_mat), na.rm = TRUE)

cellind_16wk_zscore_bycell_heatmap <- pheatmap(
  cellind_zscore_cell_heatmap_mat,
  color        = colorRampPalette(c("blue", "white", "red"))(100),
  breaks       = seq(-cellind_zscore_cell_max_abs, cellind_zscore_cell_max_abs, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  angle_col = 45,
  fontsize_row = 10,
  fontsize_col = 10,
  main         = "Metabolic pathway z-scores by cell type at 16 weeks"
)



########################################################################################
## cell independent zscore by age and cell type ----

# cell-indepedent heatmap matrix
cellind_zscore_agecell_heatmap_mat <- zscore_metadata %>%
  select(immune_vasc_ann, Age, all_of(celltype_zscore_cols)) %>%
  group_by(immune_vasc_ann, Age) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  unite("group", immune_vasc_ann, Age, sep = "_") %>%  # combine into one label
  column_to_rownames("group") %>%
  t() %>%
  as.matrix()

# 3. Clean up row names
rownames(cellind_zscore_agecell_heatmap_mat) <- gsub("_celltype_zscore", "", rownames(cellind_zscore_agecell_heatmap_mat))

# 3. Annotation bar — split on known age suffixes FOR AGE AND CELL TYPE
zscore_annotation_col <- data.frame(
  #Cell_Type = gsub("_(16|36|59|92)$", "", colnames(zscore_agecell_heatmap_mat)),
  Age_Group = gsub(".*_(16|36|59|92)$", "\\1", colnames(cellind_zscore_agecell_heatmap_mat)),
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

# save figs
ggsave(filename = paste0(date, "immunevasc_cellind_zscore_16wksbycell_heatmap_noTcells.png"), plot = cellind_16wk_zscore_bycell_heatmap, path = output_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "immunevasc_cellind_zscore_byagecell_heatmap_noTcells.png"), plot = cellind_zscore_byagecell_heatmap, path = output_path, height = 8, width = 11, dpi = 600)

ggsave(filename = paste0(date, "cellindependent_zscore_byage_heatmap.png"), plot = cellind_zscore_byage_heatmap, path = fig_path, height = 8, width = 11, dpi = 600)
