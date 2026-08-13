########################################################################
# Name: CCL_ST_reconstruct_rctdpy_WTAgingXMet_04.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 04. converting rctd-py results to rctd spacexr S4 object and RCTD post-processing
# Input: rctd-py output directory
# Final Output: post processed RCTD spacexr S4 object 
# Date created: 8/5/26
# Last updated: 8/7/26
# Author: Chloe Lucido
########################################################################

# Load Libraries ----
library(spacexr)
library(Seurat)
library(SeuratDisk)
library(future)
library(ggplot2)
library(arrow) # read in parquet 
library(hdf5r)
library(presto)
library(glmGamPoi)
library(readr)
library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)
library(patchwork)
library(sceasy)
library(reticulate)
library(SPLIT)
library(qs2) 
library(RColorBrewer)
library(Polychrome)
library(purrr)
options(future.globals.maxSize = 100 *1024^3)

# set seed
set.seed(42)

######### PATHS ##########
rctdpy_out_dir <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/RCTD_SPLIT/20260805_rctd_results_xe_Subcluster/"
output_path <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/RCTD_SPLIT/"
date <- "20260805_"

# check parquet files first 
cell_ids_df     <- read_parquet(file.path(rctdpy_out_dir, "cell_ids.parquet"))
weights_df      <- read_parquet(file.path(rctdpy_out_dir, "weights.parquet"))
spot_results_df <- read_parquet(file.path(rctdpy_out_dir, "spot_results.parquet"))
pixel_mask_df   <- read_parquet(file.path(rctdpy_out_dir, "pixel_mask.parquet"))

nrow(cell_ids_df)
nrow(weights_df)
nrow(spot_results_df)
nrow(pixel_mask_df)
sum(pixel_mask_df$pixel_mask)


# reconstruct class_df for reconstruct_rctd_from_rctdpy()
subcluster_to_broadcluster <- c(
  # "subcluster" = "broad_cluster"
  "NK Cells" = "NK Cells",
  "T Cells" = "T Cells",
  "Microglia" = "Microglia",
  "BAM" = "BAM",
  "Monocytes" = "Monocytes",
  "DC" = "DC",
  "ABC" = "ABC",
  "VLMC" = "VLMC",
  "Pericytes" = "Pericytes",
  "SMC" = "SMC",
  "Endothelial Cells" = "Endothelial Cells",
  "OPC" = "OPC",
  "Oligodendrocytes" = "Oligodendrocytes",
  "Astrocytes" = "Astrocytes",
  "Astroependymal" = "Astroependymal",
  "Tanycytes" = "Tanycytes",
  "Ependymal" = "Ependymal",
  "Hypendymal" = "Hypendymal",
  "CP" = "CP",
  "Dopaminergic_Neurons" = "Dopaminergic_Neurons",
  "Serotonergic_Neurons" = "Serotonergic_Neurons",
  "LAMP5" = "GABAergic Neuron",
  "PVALB" = "GABAergic Neuron",
  "SNCG" = "GABAergic Neuron",
  "SST" = "GABAergic Neuron",
  "VIP" = "GABAergic Neuron",
  "L2_3_IT" = "Glutamatergic Neuron",
  "L5_IT" = "Glutamatergic Neuron",
  "L6_IT" = "Glutamatergic Neuron",
  "L5_ET" = "Glutamatergic Neuron",
  "L6_CT" = "Glutamatergic Neuron"
)

# convert class_df into dataframe and set rownames
class_df <- data.frame(class = subcluster_to_broadcluster, row.names = names(subcluster_to_broadcluster))

# constructing rctd-py into a spacexr S4 object and runs SPLIT post-processing pipeline 
rctd <- reconstruct_rctd_from_rctdpy(rctdpy_out_dir, class_df = class_df)

# save rctd obj
qs_save(rctd, paste0(output_path, date, "rctd.qs2"))

