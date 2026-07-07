########################################################################
# Name: CCL_ST_SpaNorm_RPCA_WTAgingXMet_08.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 08. stats on cell proportions and cell props line graphs  
# Input Files: 
# Final Output Files: 
# Date created: 6/16/26
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
library(speckle) # for propeller


# ONLY TO BE RUN ON MCC 
options(future.globals.maxSize = 100 * 1024^3)  # slightly under your --mem

# set seed
set.seed(42)


######### Paths ############
input_obj <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/SpaNorm_RPCA/propeller_cellprops/20260616_finecelltype_props.qs2"
metadata_xlsx <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/metadata/WT_AgingXMetabolism_Xenium_metadata.xlsx"
output_path <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/"
date <- "20260619_"
###########################

# fine celltype col palette 
## fine_celltype palette ----
fine_col_palette <- c(
  # Cortical excitatory — blue gradient flowing left→right across UMAP
  "L2/3 IT neuron" = "#BBDEFB",  # leftmost, lightest sky blue
  "L2/3 IT/PVALB neuron" = "#1F78B4",  # left-center, strong blue
  "L5 IT neuron" = "#7986CB",  # center, blue-indigo
  "L5/6 NP neuron" = "#5C6BC0",  # center, deeper indigo
  "L6 CT neuron" = "#CAB2D6",  # center-right, soft lavender
  "L6 CT/PVALB neuron" = "#7E57C2",  # right, deep purple
 
  # Inhibitory — warm pinks/mauves, overlapping cortical blob
  "VIP neuron" = "#F9C6D0",
  "LAMP5 neuron" = "#D81B60",
  "SST neuron" = "#CE93D8",

  # Detached neurons — spatially isolated, distinct hues
  "Thalamic neuron" = "#33A02C",  # ← was #D4A800
    "Dentate Gyrus neuron" = "#D4A800",  # ← was #33A02C
 "Hippocampal neuron" = "#BCCE6B",  # upper-right near DG

  # Vascular — oranges (unchanged)
  "Endothelial cell" = "#FDBF6F",
  "Vascular mural" = "#FF7F00",
  "Pericyte" = "#C45E00",

  # Glia & other (unchanged)
  "Astrocyte" = "#E64B35",
  "Oligodendrocyte" = "#B2DF8A",
  "OPC" = "#A8E6CF",
  "Microglia" = "#00897B",
  "Ependymal cell" = "#FB9A99",
  "Choroid plexus" = "#5C8A3C",
  "T cell" = "#8B5E3C"
)

# category map for line graphs ----
# ---- Step 1: Add category labels to summary_df ----
category_map <- data.frame(
  CellType = c(
    # exc cortical neuronal types 
    "L2/3 IT neuron", "L2/3 IT/PVALB neuron", "L5 IT neuron", "L5/6 NP neuron", "L6 CT neuron", "L6 CT/PVALB neuron", 
    
    # inh cortical neu
    "VIP neuron", "LAMP5 neuron", "SST neuron", 
    
    # other neu
    "Thalamic neuron", "Dentate Gyrus neuron", "Hippocampal neuron",
    
    # immune
    "Microglia", "T cell", 
    
    # vascular-related
    "Endothelial cell", "Pericyte", "Vascular mural", "Ependymal cell", 
    #"Choroid plexus", # removing CP 
               
    # other glia
    "Astrocyte", "Oligodendrocyte", "OPC"
                 ),  
  
  
  Category = c("Neuronal", "Neuronal", "Neuronal", "Neuronal", "Neuronal", "Neuronal", 
               "Neuronal", "Neuronal", "Neuronal", 
               "Neuronal", "Neuronal", "Neuronal",
               "Immune", "Immune", 
               "Vascular-related", "Vascular-related", "Vascular-related", "Vascular-related", 
               #"Vascular-related", 
               "Glia", "Glia", "Glia")  
)



# CELL PROPORTIONS WITH PROPELLER----

fine_props <- qs_read(input_obj)
metadata <- read_xlsx(metadata_xlsx)



# The raw (untransformed) proportions are stored here:
raw_fine_props <- fine_props$Proportions  # rows = cell types, columns = samples

# ---- Step 2: Add Age group info per sample ----
# Build a sample -> age lookup from your Seurat metadata
sample_age <- metadata %>%
  select(sample_ID, Age) %>%
  distinct()

# Reshape proportions to long format
fine_props_long <- as.data.frame(raw_fine_props) %>%
  rename(CellType = clusters, sample_ID = sample, Proportion = Freq)

# create age column using sample IDs
fine_props_long <- fine_props_long %>%
  mutate(Age = as.numeric(gsub("wk.*", "", sample_ID)))  # extracts 16, 36, 59, 92


# ---- Step 3: Calculate Mean, SD, SE per cell type per age ----
fine_summary_df <- fine_props_long %>%
  group_by(CellType, Age) %>%
  summarise(
    Mean = mean(Proportion),
    SD   = sd(Proportion),
    N    = n(),
    SE   = SD / sqrt(N),
    .groups = "drop"
  )

# add category map to fine_summary_df 
fine_summary_df <- fine_summary_df %>%
  left_join(category_map, by = "CellType")

print(fine_summary_df)


# ---- Significance labels from your propeller output ----
# Add all your cell types and their FDR values here
fine_sig_df <- data.frame(
  CellType = c("Microglia", "Oligodendrocyte", "T cell", "Endothelial cell"),        # MANUALLY ADD
  FDR      = c(0.001002369, 0.003747815, 0.003747815, 0.026942719)    # MANUALLY ADD
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


plot_category_norm <- function(data, category_name) {
  plot_data <- data %>% filter(Category == category_name)
  category_colors <- fine_col_palette[names(fine_col_palette) %in% unique(plot_data$CellType)]
  category_colors <- category_colors[order(names(category_colors))]
  
  line_end_labels <- plot_data %>%
    filter(Age == 92, !is.na(sig_label), sig_label != "ns")
  
  ggplot(plot_data, aes(x = Age, y = Mean_norm, color = CellType, fill = CellType, group = CellType)) +
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
    scale_color_manual(values = category_colors) +
    scale_fill_manual(values = category_colors) +
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
                            fill = category_colors,
                            color = "white", 
                            stroke = 0.3,
                            linetype = 0)
      ),
      fill = "none"    # suppress the fill legend
    ) +
    labs(
      title = category_name,
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
p_neuronal <- plot_category_norm(fine_summary_df_normalized, "Neuronal")
p_glia     <- plot_category_norm(fine_summary_df_normalized, "Glia")
p_vascular <- plot_category_norm(fine_summary_df_normalized, "Vascular-related")
p_immune   <- plot_category_norm(fine_summary_df_normalized, "Immune")

# View individually
p_neuronal # saving with geom_ribbon alpha = 0.1
p_glia
p_vascular
p_immune

ggsave(paste0(output_path, date, "norm_to_16wks_neuronal_props_linegraph.png"), plot = p_neuronal, width = 11, height = 8, dpi = 600)
ggsave(paste0(output_path, date, "norm_to_16wks_glia_props_linegraph.png"), plot = p_glia, width = 11, height = 8, dpi = 600)
ggsave(paste0(output_path, date, "norm_to_16wks_vasc_NO_CP_props_linegraph.png"), plot = p_vascular, width = 11, height = 8, dpi = 600)
ggsave(paste0(output_path, date, "norm_to_16wks_immune_props_linegraph.png"), plot = p_immune, width = 11, height = 8, dpi = 600)
