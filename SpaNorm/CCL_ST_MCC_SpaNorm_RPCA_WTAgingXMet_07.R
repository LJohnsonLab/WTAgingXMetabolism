########################################################################
# Name: CCL_ST_MCC_SpaNorm_RPCA_WTAgingXMet_07.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 06. making UMAPs, imgdimplots, and calculating cell proportions btwn groups for both gen and fine celltypes  
# Input Files: 
# Final Output Files: 
# Date created: 6/14/26
# Last updated: 7/7/26
# Author: Chloe Lucido
########################################################################

# Auto-install missing packages
required_packages <- c("Seurat", "future", "ggplot2", "readr", "dplyr", 
                       "data.table","tidyr", "tibble", "stringr", "forcats", "lubridate",  "readxl", "patchwork", 
                       "qs2", "RColorBrewer", "Polychrome", "purrr",
                       "presto", "glmGamPoi", "SCP", "pheatmap")

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
library(speckle) # for propeller (cell prop stats)

# ONLY TO BE RUN ON MCC 

# Use the cores you requested
plan("multicore", workers = 8)
options(future.globals.maxSize = 400 * 1024^3)  # slightly under your --mem

# set seed
set.seed(42)

# to prevent this error: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so.
options(future.rng.onMisuse = "resolve")
options(future.seed = TRUE)


######### Paths ############
input_obj <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260612_SpaNorm_RPCA_cell_region_ann.qs2"
fig_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/SpaNorm_RPCA_figs/"
output_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/"
cellprops_outputpath <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/propeller_cellprops/"
date <- "20260623_"
###########################


# 01. Read in obj ----
SpaNorm_RPCA.obj <- qs_read(input_obj)

# fix L3 ET neuron name
Idents(SpaNorm_RPCA.obj) <- "fine_celltype"

# NOTE:
## after looking at  split image dim plot, i realized that what I annotated as L3 ET neurons is actually hippocampal neurons
## changing this annotation to reflect this
SpaNorm_RPCA.obj <- RenameIdents(SpaNorm_RPCA.obj, 
                                 "L3 ET Neuron" = "Hippocampal neuron")

SpaNorm_RPCA.obj$fine_celltype <- Idents(SpaNorm_RPCA.obj)


# fix cluster order for aesthetic reasons ----
## reorder gen_cell type ----
SpaNorm_RPCA.obj$gen_celltype <- factor(SpaNorm_RPCA.obj$gen_celltype, 
                                        levels = c(
                                                   #neurons
                                                   "Neuron", 
                                                   # glia
                                                    "Astrocyte",
                                                    "Oligodendrocyte", 
                                                    "OPC", 
                                                    # immune
                                                    "Microglia", 
                                                    "T cell", 
                                                    # vascular
                                                    "Endothelial cell", 
                                                    "Pericyte",
                                                    "Vascular mural", 
                                                    "Ependymal cell",
                                                    "Choroid plexus")) 

## reorder fine_cell type ----

SpaNorm_RPCA.obj$fine_celltype <- factor(SpaNorm_RPCA.obj$fine_celltype, 
                                        levels = c(
                                                    # neurons
                                                    "L2/3 IT neuron",
                                                    "L2/3 IT/PVALB neuron",
                                                    "L5 IT neuron",
                                                    "L5/6 NP neuron",
                                                    "L6 CT neuron",
                                                    "L6 CT/PVALB neuron",
                                                    "LAMP5 neuron",
                                                    "SST neuron",
                                                    "VIP neuron",

                                                    "Thalamic neuron",
                                                    "Dentate Gyrus neuron",
                                                    "Hippocampal neuron",
                                                   # glia
                                                    "Astrocyte",
                                                    "Oligodendrocyte", 
                                                    "OPC", 
                                                    # immune
                                                    "Microglia", 
                                                    "T cell", 
                                                    # vascular
                                                    "Endothelial cell", 
                                                    "Pericyte",
                                                    "Vascular mural", 
                                                    "Ependymal cell",
                                                    "Choroid plexus")) 


# cluster palettes ----
## gen_celltype palette ----
gen_col_palette <- c(
                     "Neuron" = "#1F78B4",
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

# SANITY CHECK
#levels(SpaNorm_RPCA.obj$gen_celltype)

# [1] "Neuron"           "Astrocyte"        "Oligodendrocyte"  "OPC"             
# [5] "Microglia"        "T cell"           "Endothelial cell" "Pericyte"        
# [9] "Vascular mural"   "Ependymal cell"   "Choroid plexus"  

#table(SpaNorm_RPCA.obj$gen_celltype)
#          Neuron        Astrocyte  Oligodendrocyte              OPC 
#          739646           194195           304879            47206 
#       Microglia           T cell Endothelial cell         Pericyte 
#           53886             1160           106701            29778 
#  Vascular mural   Ependymal cell   Choroid plexus 
#           94088            21532            12009 

  
library(SCP) # for pretty UMAPs
# make figures 

## gen cell type UMAP ---- 
gen_UMAP <- CellDimPlot(
  srt = SpaNorm_RPCA.obj,
  group.by = "gen_celltype",
  reduction = "umap",
  theme_use = "theme_blank",
  pt.size = 0.6
) + scale_color_manual(values = gen_col_palette) +
                      guides(color = guide_legend(title = NULL, # remove legend title
                                                    ncol = 1, # keep legend in 1 column
                                                    override.aes = list(size = 3))) # legend dot size

## fine_celltype UMAP ----
fine_UMAP <- CellDimPlot(
  srt = SpaNorm_RPCA.obj,
  group.by = "fine_celltype",
  reduction = "umap",
  theme_use = "theme_blank",
  pt.size = 0.6
) + scale_color_manual(values = fine_col_palette) + 
                        guides(color = guide_legend(title = NULL, # remove legend title
                                                    ncol = 1, # keep legend in 1 column
                                                    override.aes = list(size = 3))) # legend dot size

## gen grouped img dim plot ----
gen_imagedim_fov1 <- ImageDimPlot(SpaNorm_RPCA.obj,
                              group.by = "gen_celltype",
                              cols = gen_col_palette,
                              dark.background = F,
                              fov = "fov") + labs(fill = NULL)

## gen split imgdimplot ----
gen_imagedim_splitfov1 <- ImageDimPlot(SpaNorm_RPCA.obj,
                              group.by = "gen_celltype", 
                                split.by = "gen_celltype",
                              cols = gen_col_palette,
                              dark.background = F,
                              fov = "fov") + labs(fill = NULL)


## fine grouped imgdimplot ----
fine_imagedim_fov1 <- ImageDimPlot(SpaNorm_RPCA.obj,
                              group.by = "fine_celltype",
                              cols = fine_col_palette,
                              dark.background = F,
                              fov = "fov") + labs(fill = NULL)

# look at images/image names that the object has 
#Images(SpaNorm_RPCA.obj)

# [1] "fov"      "fov.2"    "fov.2.3"  "fov.2.4"  "fov.2.5"  "fov.3"
# [7] "fov.3.7"  "fov.3.8"  "fov.4"    "fov.4.10" "fov.4.11" "fov.5"
# [13] "fov.5.13" "fov.5.14" "fov.5.15" "fov.6"    "fov.6.17" "fov.6.18"
# [19] "fov.7"    "fov.7.20" "fov.7.21" "fov.7.22" "fov.7.23" "fov.7.24"

#table(SpaNorm_RPCA.obj$sample_ID, SpaNorm_RPCA.obj$slide_ID)
 
#         0021991 0021998 0021999 0022573 0069039 0069045 0069118
#16wk_F1       0   66873       0       0       0       0       0
#16wk_F2       0       0       0       0       0       0   64244
#16wk_F4       0       0       0       0       0   63103       0
#16wk_M1   69680       0       0       0       0       0       0
#16wk_M2       0   71756       0       0       0       0       0
#16wk_M4       0       0       0       0   58033       0       0
#36wk_F1   67485       0       0       0       0       0       0
#36wk_F3       0       0       0       0       0   67618       0
#36wk_F4       0       0       0       0       0       0   67001
#36wk_M1       0   65784       0       0       0       0       0
#36wk_M2   75848       0       0       0       0       0       0
#36wk_M4       0       0       0       0       0       0   57099
#59wk_F2       0       0       0   80830       0       0       0
#59wk_F3       0       0       0       0       0       0   68186
#59wk_F5       0       0       0       0   49701       0       0
#59wk_M1       0       0       0   64001       0       0       0
#59wk_M3       0       0       0       0       0   64106       0
#59wk_M5       0       0       0       0       0       0   71959
#92wk_F2       0       0       0   72664       0       0       0
#92wk_F3       0       0       0       0   66941       0       0
#92wk_F5       0       0       0       0       0       0   70278
#92wk_M1       0       0   54198       0       0       0       0
#92wk_M2       0       0       0   73868       0       0       0
#92wk_M3       0       0       0       0   73824       0       0



# which FOV belongs to which samples? 
# List all FOVs and which sample_IDs they contain
for (img_name in Images(SpaNorm_RPCA.obj)) {
  cells_in_fov <- Cells(SpaNorm_RPCA.obj@images[[img_name]])
  sample_ids   <- unique(SpaNorm_RPCA.obj$sample_ID[cells_in_fov])
  cat(img_name, "→", sample_ids, "\n")
}

#fov → 92wk_M1
#fov.2 → 59wk_F2
#fov.2.3 → 92wk_M2
#fov.2.4 → 59wk_M1
#fov.2.5 → 92wk_F2
#fov.3 → 36wk_M2
#fov.3.7 → 36wk_F1
#fov.3.8 → 16wk_M1
#fov.4 → 36wk_M1
#fov.4.10 → 16wk_M2
#fov.4.11 → 16wk_F1
#fov.5 → 92wk_M3
#fov.5.13 → 92wk_F3
#fov.5.14 → 59wk_F5
#fov.5.15 → 16wk_M4
#fov.6 → 16wk_F4
#fov.6.17 → 59wk_M3
#fov.6.18 → 36wk_F3
#fov.7 → 92wk_F5
#fov.7.20 → 36wk_F4
#fov.7.21 → 59wk_M5
#fov.7.22 → 59wk_F3
#fov.7.23 → 36wk_M4
#fov.7.24 → 16wk_F2


# imagedimplot each fov to find 16wk_F2
fov7_imgdim <- ImageDimPlot(SpaNorm_RPCA.obj, 
             group.by = "fine_celltype", 
             cols = fine_col_palette, 
             dark.background = F, 
             fov = "fov.7") +
                      labs(fill = NULL)

fov7_20_imgdim <- ImageDimPlot(SpaNorm_RPCA.obj,
              group.by = "fine_celltype",
              cols = fine_col_palette,
              dark.background = F,
              fov = "fov.7.20") +
                       labs(fill = NULL)


fov7_21_imgdim <- ImageDimPlot(SpaNorm_RPCA.obj,
              group.by = "fine_celltype",
              cols = fine_col_palette,
              dark.background = F,
              fov = "fov.7.21") +
                       labs(fill = NULL)



fov7_22_imgdim <- ImageDimPlot(SpaNorm_RPCA.obj,
              group.by = "fine_celltype",
              cols = fine_col_palette,
              dark.background = F,
              fov = "fov.7.22") +
                       labs(fill = NULL)

fov7_23_imgdim <- ImageDimPlot(SpaNorm_RPCA.obj,
              group.by = "fine_celltype",
              cols = fine_col_palette,
              dark.background = F,
              fov = "fov.7.23") +
                       labs(fill = NULL)

                   
fov7_24_imgdim <- ImageDimPlot(SpaNorm_RPCA.obj,
              group.by = "fine_celltype",
              cols = fine_col_palette,
              dark.background = F,
              fov = "fov.7.24") +
                       labs(fill = NULL)


fov7_24_imgdim_split <- ImageDimPlot(SpaNorm_RPCA.obj,
              group.by = "fine_celltype",
              split.by = "fine_celltype",
              cols = fine_col_palette,
              dark.background = F,
              fov = "fov.7.24") +
                       labs(fill = NULL)

# save figs 
ggsave(filename = paste0(date, "fov7_imgdim.png"), plot = fov7_imgdim, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "fov7_20_imgdim.png"), plot = fov7_20_imgdim, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "fov7_21_imgdim.png"), plot = fov7_21_imgdim, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "fov7_22_imgdim.png"), plot = fov7_22_imgdim, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "fov7_23_imgdim.png"), plot = fov7_23_imgdim, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "fov7_24_imgdim.png"), plot = fov7_24_imgdim, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "fov7_24_imgdim_split.png"), plot = fov7_24_imgdim_split, path = fig_path, height = 8, width = 11, dpi = 600)

# subsetting to cells in rep brain (16wk_F2) 
rep.obj <- subset(SpaNorm_RPCA.obj, subset = sample_ID == "16wk_F2")

# reattaching image
rep.obj@images[["fov.7"]] <- SpaNorm_RPCA.obj@images[["fov.7"]]

# Then subset the FOV to only the cells in rep.obj
rep.obj@images[["fov.7"]] <- subset(rep.obj@images[["fov.7"]], 
#                                     cells = Cells(rep.obj))
# Then verify images transferred
Images(rep.obj)

## rep fine imagedim plot of 16wk F2 brain ----
rep_imagedim <- ImageDimPlot(rep.obj,
                              group.by = "fine_celltype",
                              cols = fine_col_palette,
                              dark.background = F
                             # ,fov = "fov.7"
                              ) + labs(fill = NULL)


## gen split imgdimplot ----
fine_imagedim_splitfov1 <- ImageDimPlot(SpaNorm_RPCA.obj,
                              group.by = "fine_celltype",
                                        split.by = "fine_celltype",
                              cols = fine_col_palette,
                              dark.background = F,
                              fov = "fov") + labs(fill = NULL)


# dimplots
ggsave(filename = paste0(date, "gen_cell_UMAP.png"), plot = gen_UMAP, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "fine_cell_UMAP.png"), plot = fine_UMAP, path = fig_path, height = 8, width = 11, dpi = 600)

ggsave(filename = paste0(date, "gen_imagedim_fov1.png"), plot = gen_imagedim_fov1, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "fine_imagedim_fov1.png"), plot = fine_imagedim_fov1, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "gen_imagedim_splitfov1.png"), plot = gen_imagedim_splitfov1, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "fine_imagedim_splitfov1.png"), plot = fine_imagedim_splitfov1, path = fig_path, height = 8, width = 11, dpi = 600)
ggsave(filename = paste0(date, "rep_imagedim.png"), plot = rep_imagedim, path = fig_path, height = 8, width = 11, dpi = 600)


# save obj with updated factor levels and updated L3 ET name
qs_save(SpaNorm_RPCA.obj, paste0(output_path, date, "SpaNorm_RPCA.qs2"))

# Propeller: cell proportion stats btwn age groups ----
## gen propeller stats ----
gen_propeller_results <- propeller(
  x        = SpaNorm_RPCA.obj,
  clusters = SpaNorm_RPCA.obj$gen_celltype,
  group    = SpaNorm_RPCA.obj$Age,
  sample   = SpaNorm_RPCA.obj$sample_ID,
  trend    = FALSE,
  robust   = TRUE,
  transform = "logit"
)

## fine propeller stats ----
fine_propeller_results <- propeller(
  x        = SpaNorm_RPCA.obj,
  clusters = SpaNorm_RPCA.obj$fine_celltype,
  group    = SpaNorm_RPCA.obj$Age,
  sample   = SpaNorm_RPCA.obj$sample_ID,
  trend    = FALSE,
  robust   = TRUE,
  transform = "logit"
)


# create celltype column before saving as csv 
gen_propeller_results <- gen_propeller_results %>% 
    rownames_to_column("CellType")

fine_propeller_results <- fine_propeller_results %>%
    rownames_to_column("CellType")

# save propeller stats ----
write_csv(gen_propeller_results, paste0(cellprops_outputpath, date, "gencelltype_props_stats.csv"))
write_csv(fine_propeller_results, paste0(cellprops_outputpath, date, "finecelltype_props_stats.csv"))



gen_props <- getTransformedProps(
  clusters = SpaNorm_RPCA.obj$gen_celltype,
  sample   = SpaNorm_RPCA.obj$sample_ID,
  transform = "logit"   # match what you used in propeller
)

fine_props <- getTransformedProps(
  clusters = SpaNorm_RPCA.obj$fine_celltype,
  sample   = SpaNorm_RPCA.obj$sample_ID,
  transform = "logit"   # match what you used in propeller
)

# save props objs as qs2 files 
qs_save(gen_props, paste0(cellprops_outputpath, date, "gencelltype_props.qs2"))

qs_save(fine_props, paste0(cellprops_outputpath, date, "finecelltype_props.qs2"))
