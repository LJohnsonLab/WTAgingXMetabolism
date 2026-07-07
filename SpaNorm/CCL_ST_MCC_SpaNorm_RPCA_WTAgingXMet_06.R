########################################################################
# Name: CCL_ST_SpaNorm_RPCA_WTAgingXMet_06.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - SpaNorm RPCA 06. renaming clusters to cell types, and adding regional annotations 
# Input Files:   
# Final Output Files: 
# Date created: 6/11/26
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
input_obj <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260609_SpaNorm_RPCAint_processed.qs2" # not using obj with subclusters metadata bc I don't need that 
output_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/"
fig_output_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/SpaNorm_RPCA_figs/"
XE_ann_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/XE_regional_annotations/"
date <- "20260612_"

#######################
# NOTES: 
## first half of script (gen_celltype and regional annotations addns) done on 6/12/26
## second half of script (fine_celltype annotations) done on 6/15/26
## did not rename the final obj because it did not majorly change 

# 01. Read in SpaNorm RPCA processed obj (no subcluster metadata bc i am not using them) ---- 
SpaNorm_RPCA.obj <- qs_read(input_obj)

# 02. Annotating Cell types ---- 
##02a. make sure idents are correct ----
Idents(SpaNorm_RPCA.obj) <- "clusters_res05"

## 02b. Rename Idents to GENERAL cell types ----
SpaNorm_RPCA.obj <- RenameIdents(SpaNorm_RPCA.obj, 
                                 "0" = "Neuron", 
                                 "1" = "Oligodendrocyte", 
                                 "2" = "Endothelial cell", 
                                 "3" = "Astrocyte", 
                                 "4" = "Oligodendrocyte", 
                                 "5" = "Neuron", 
                                 "6" = "Astrocyte", 
                                 "7" = "Neuron", 
                                 "8" = "Neuron", 
                                 "9" = "Vascular mural", 
                                 "10" = "Oligodendrocyte", 
                                 "11" = "Neuron", 
                                 "12" = "Microglia", 
                                 "13" = "Neuron",
                                 "14" = "Neuron", 
                                 "15" = "OPC", 
                                 "16" = "Neuron", 
                                 "17" = "Neuron", 
                                 "18" = "Oligodendrocyte", 
                                 "19" = "Vascular mural", 
                                 "20" = "Neuron",
                                 "21" = "Pericyte", 
                                 "22" = "Neuron", 
                                 "23" = "Ependymal cell", 
                                 "24" = "Neuron", 
                                 "25" = "Neuron", 
                                 "26" = "Choroid plexus", 
                                 "27" = "Neuron", 
                                 "28" = "Neuron", 
                                 "29" = "Neuron", 
                                 "30" = "TRASH", 
                                 "31" = "T cell", 
                                 "32" = "TRASH", 
                                 "33" = "TRASH")


## 02c. subset OUT trash clusters ----
SpaNorm_RPCA.obj <- subset(SpaNorm_RPCA.obj, idents = "TRASH", invert = T)

# Note:
## combined, clusters 30, 32 and 33 have 4328 cells
## these clusters make up ~0.2% of all cells so I am not renormalizing or scaling

## 02d. quick dimplot for annotation confirmation -----
UMAP_1 <- DimPlot(SpaNorm_RPCA.obj, label = T, cols = "polychrome")

ggsave(filename = paste0(date, "SpaNorm_RPCA_gen_celltype_UMAP.png"), plot = UMAP_1, path = fig_output_path, width = 11, height = 10)


## 02e. save gen celltype annotations to metadata ----
SpaNorm_RPCA.obj$gen_celltype <- Idents(SpaNorm_RPCA.obj)


# 03. add regional annotations to metadata ----

## 03a. read in XE regional annotation files ----

files <- list.files(
  path = XE_ann_path, 
  pattern = "cells_stats\\.csv$", # grab files containing "cell_stats" in file name and the extension .csv
  recursive = T, # automatically searches subfolders 
  full.names = T # returns complete file paths 
 

## 03b. create annotation list ---- 
annotation_list <- map(files, function(path) {

  # read file
  df <- read.csv(path, skip = 2, blank.lines.skip = TRUE) # first 2 lines do not contain necessary information
  
  # extract filename
  fname <- basename(path)
  
  # Extract slide number (9 digits)
  slide_number <- str_extract(fname, "(?<=Run\\d_)\\d+")
  
  # Extract sample ID (age + sex, e.g. "30wk_M9")
  sample_id <- str_extract(fname, "\\d+wk_[MF]\\d+")
  
  # Extract region prefix (e.g. "HYP")
  region <- str_extract(fname, "(?<=_)[A-Z]+(?=_cells_stats)")
  
  # Determine slide prefix from last 3 digits of slide number
  slide_suffix <- str_sub(slide_number, -3, -1)
  slide_prefix <- paste0("slide", slide_suffix, "_")
  
  # add these metadata columns to the df containing cell stats annotation info
  df %>%
    mutate(
      Cell.ID = paste0(slide_prefix, Cell.ID),
      sample_ID = sample_id,
      anatomical_region = region, 
      slide = slide_number
    )
})

# print annotation_list colnames 
print(colnames(annotation_list[[1]]))

# check if any files are causing issues 
map(annotation_list, dim)

# check colnames of each element
map(annotation_list, colnames)

## 03c. format information to add to seurat obj metadata ----

# bind all rows from annotation_list into one big df
annotation_df <- bind_rows(annotation_list)

# make sure to keep only unique cells (gets rid of duplicating cells)
annotation_df_unique <- annotation_df %>%
  distinct(Cell.ID, .keep_all = TRUE)

# set rownames to cell IDs to add to seurat metadata 
rownames(annotation_df_unique) <- annotation_df_unique$Cell.ID

annotation_df_unique <- annotation_df_unique[Cells(SpaNorm_RPCA.obj), , drop = FALSE] # CHANGE OBJ

## 03d. add to seurat obj metadata ----
SpaNorm_RPCA.obj <- AddMetaData(SpaNorm_RPCA.obj, metadata = annotation_df_unique["anatomical_region"]) # changed it from "Region" bc of squidpy analyses



# 04. save obj w/ gen_celltype and regional annotation metadata cols ----
qs_save(SpaNorm_RPCA.obj, paste0(output_path, date, "SpaNorm_RPCA_cell_region_ann.qs2"))



# 6/15/2026 adding fine_celltype metadata col for ICBEM ##############################
# 01. read in file----
SpaNorm_RPCA.obj<- qs_read(paste0(output_path, date, "SpaNorm_RPCA_cell_region_ann.qs2"))


# 02. add fine_celltype annotations (neuronal types) ----
## 02a.  make sure idents are correct ----
Idents(SpaNorm_RPCA.obj) <- "clusters_res05"

print(unique(SpaNorm_RPCA.obj$clusters_res05))

## 02b. Rename Idents to FINE cell types ----
SpaNorm_RPCA.obj <- RenameIdents(SpaNorm_RPCA.obj,
                                  "0" = "VIP neuron",
                                  "1" = "Oligodendrocyte",
                                  "2" = "Endothelial cell",
                                  "3" = "Astrocyte",
                                  "4" = "Oligodendrocyte",
                                  "5" = "LAMP5 neuron",
                                  "6" = "Astrocyte",
                                  "7" = "L6 CT neuron",
                                  "8" = "Thalamic neuron",
                                  "9" = "Vascular mural",
                                  "10" = "Oligodendrocyte",
                                  "11" = "L2/3 IT neuron",
                                  "12" = "Microglia",
                                  "13" = "Dentate Gyrus neuron",
                                  "14" = "L2/3 IT/PVALB neuron",
                                  "15" = "OPC",
                                  "16" = "L5 IT neuron",
                                  "17" = "SST neuron",
                                  "18" = "Oligodendrocyte",
                                  "19" = "Vascular mural",
                                  "20" = "L5/6 NP neuron",
                                  "21" = "Pericyte",
                                  "22" = "L3 ET neuron",
                                  "23" = "Ependymal cell",
                                  "24" = "L5/6 NP neuron",
                                  "25" = "Hippocampal neuron",
                                  "26" = "Choroid plexus",
                                  "27" = "LAMP5 neuron",
                                  "28" = "VIP neuron",
                                  "29" = "L6 CT/PVALB neuron",
                                  "31" = "T cell"
                                  )
 
 
# dont need to subset OUT trash clusters bc i already did that in this obj


## 02c. add fine_celltype metadata col ----
SpaNorm_RPCA.obj$fine_celltype <- Idents(SpaNorm_RPCA.obj)

## 02d. quick dimplot for annotation confirmation -----
UMAP_2 <- DimPlot(SpaNorm_RPCA.obj, label = T, cols = "polychrome")

ggsave(filename = paste0(date, "SpaNorm_RPCA_fine_celltype_UMAP.png"), plot = UMAP_2, path = fig_output_path, width = 11, height = 10)

# 03. save updated file ----
qs_save(SpaNorm_RPCA.obj, paste0(output_path, date, "SpaNorm_RPCA_cell_region_ann.qs2"))


# making sure all the metadata cols that I added are good

table(SpaNorm_RPCA.obj$clusters_res05)

#     0      1      2      3      4      5      6      7      8      9     10 
#163584 116327 106701 102890  93115  91456  91305  63755  62690  62497  59685 
#    11     12     13     14     15     16     17     18     19     20     21 
# 59006  53886  48118  48022  47206  44959  40801  35752  31591  30399  29778 
#    22     23     24     25     26     27     28     29     30     31     32 
# 28447  21532  17016  14184  12009  11819   9188   6202      0   1160      0 
#    33 
#     0


table(SpaNorm_RPCA.obj$gen_celltype)

#          Neuron  Oligodendrocyte Endothelial cell        Astrocyte 
#          739646           304879           106701           194195 
#  Vascular mural        Microglia              OPC         Pericyte 
#           94088            53886            47206            29778 
#  Ependymal cell   Choroid plexus           T cell 
#           21532            12009             1160 


table(SpaNorm_RPCA.obj$fine_celltype)


#          VIP neuron      Oligodendrocyte     Endothelial cell 
#              172772               304879               106701 
#           Astrocyte         LAMP5 neuron         L6 CT neuron 
#              194195               103275                63755 
#     Thalamic neuron       Vascular mural       L2/3 IT neuron 
#               62690                94088                59006 
#           Microglia Dentate Gyrus neuron L2/3 IT/PVALB neuron 
#               53886                48118                48022 
#                 OPC         L5 IT neuron           SST neuron 
#               47206                44959                40801 
#      L5/6 NP neuron             Pericyte         L3 ET Neuron 
#               47415                29778                28447 
#      Ependymal cell   Hippocampal neuron       Choroid plexus 
#               21532                14184                12009 
#  L6 CT/PVALB neuron               T cell 
#                6202                 1160 
