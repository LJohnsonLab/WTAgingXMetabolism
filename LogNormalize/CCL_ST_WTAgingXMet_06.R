########################################################################
# Name: CCL_ST_WTAgingXMet_06.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 06. Metabolic Pathway z scoring on all cells in log normalize obj => to generate quick figs for Trainee seminar pres on 6/2/26 
#             - 01 
#             - 02 
#             - 03 
#             - 04 
#             - 05  
#             - 06  
#             - 07  
#             - 08  
# Input Files: final annotated obj (normalize and scale) from part 3 => 20260428_mergedobj_gen_annotations_03.qs2
# Final Output Files: 
# Date created: 5/28/26
# Last updated: 6/10/26
# Author: Chloe Lucido
########################################################################

# Load Libraries ----
library(spacexr)
library(Seurat)
library(SeuratDisk)
library(future)
library(ggplot2)
library(arrow)
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
library(readxl)
library(BPCells)
options(future.globals.maxSize = 100 *1024^3)

# set seed
set.seed(42)


# 01. read in final annotated (normalize and scale) obj from part 3 ----
normscale.obj <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260428_mergedobj_gen_annotations_03.qs2")


# 02. Metabolic pathway reference tables ----

## 02a. Fuel pathways reference table

fuel_path_table <- list(
  "Glycolysis" = c("Hk1", "Hk2", "Hk3", "Hkdc1", "Pfkl", "Pfkm", "Pfkp", "Pkm", "Ldha", "Ldhb", "Pgm1", "Pdhb"), 
  "OXPHOS" = c("Ndufv2", "Ndufs1", "Sdha", "Uqcrc1", "Cox10", "Sdhc", "Sdhb", "Cox14", "Cox4l1"), 
  "TCA" = c("Cs", "Acly", "Idh2", "Dlst", "Pdha1"), 
  "PPP" = c("G6pd2", "Tkt", "Taldo1", "Pgls", "Prps2", "Pgd", "Pgm2"), 
  "PYR_PUR" = c("Pnp", "Dpyd", "Ctps2", "Polr1a", "Impdh2", "Adss", "Ppat", "Cad", "Nt5c", "Gmps", "Adcy5"), 
  "Ketone_Body" = c("Hmgcs1", "Hmgcl", "Oxct1", "Acat1"), 
  "Pyruvate" = c("Acss2", "Dlat", "Acaca", "Pcx", "Mdh1", "Mdh2", "Me1", "Aldh2"), 
  "FRUC_MAN" = c("Pmm1", "Pfkfb3", "Khk", "Aldob"), 
  "STA_SUC" = c("Gaa", "Ganc", "Pygb", "G6pc3", "Gys2", "Gck"), 
  "PANT_COA" = c("Pank3", "Coasy", "Ppcdc"), 
  "Fatty_Acid" = c("Acsl6", "Acsl4", "Hadhb", "Ppt1", "Hsd17b12", "Hacd2", "Tecr", "Cpt1b", "Fasn", "Abhd17a", "Acadl", "Ecl1", "Elovl1"), 
  "Glutathione" = c("Gclc", "Gclm", "Lap3", "Gstm1", "Idh1")
  )

## 02b. Amino Acid pathways reference table ----

AA_path_table <- list(
  "Ala_Asp_Glu" = c("Got1", "Gpt2", "Glud1", "Got2", "Folh1", "Aspa"), 
  "Gly_Ser_Thr" = c("Phgdh", "Dld", "Gatm", "Shmt2", "Pgam1", "Bpgm"), 
  "Cys_Met" = c("Bckdha", "Hibadh", "Auh", "Cbs", "Acad8"), 
  "Lys" = c("Ogdh", "Hadha", "Echs1"), 
  "Arg_Pro" = c("Odc1", "Sat1", "Prodh", "P4ha1", "Pycr2")
  #, 
  #"His" = c(), 
  #"Tyr" = c(), 
  #"Phe" = c(), 
  #"Tryp" = c()
) 

# 03. trying z-score method ---- 

library(Matrix)

fuelpathway_scores <- imap(fuel_path_table, function(genes, pathway_name) {
  
  # Step 1: get genes present in object
  pathway_genes <- intersect(genes, rownames(normscale.obj))
  
  # Step 2: extract sparse matrix (genes x cells) - no dense conversion
  mat <- GetAssayData(normscale.obj, layer = "data")[pathway_genes, , drop = F]
  
  # Step 3: z-score each gene across cells (row-wise), keeping sparse structure
  # Scale each gene: subtract mean, divide by sd
  row_means <- rowMeans(mat)
  row_sds   <- apply(mat, 1, sd)
  row_sds[row_sds == 0] <- 1  # avoid divide-by-zero for zero-variance genes

  
  
  colMeans((mat - row_means) / row_sds)  # broadcasts correctly
  
}) |> 
  compact() |> # drops nulls from skipped pathways 
  as.data.frame() |> 
  rename_with(~ paste0(.x, "_zscore")) # append zscore to pathway cols names 


# add pathway z scores for each cell to metadata 
normscale.obj <- AddMetaData(normscale.obj, metadata = fuelpathway_scores)



# Visualise
fuelpath_groupedby_age_dotplot <- DotPlot(normscale.obj, 
        features = c(
          #"Glycolysis_zscore", 
                     "OXPHOS_zscore"
                     #, 
                     #"TCA_zscore", 
                     #"PPP_zscore"
                     #, 
                     #"PYR_PUR_zscore", 
                     #"Ketone_Body_zscore", 
                     #"Pyruvate_zscore", 
                     #"FRUC_MAN_zscore", 
                     #"STA_SUC_zscore", 
                     #"PANT_COA_zscore", 
                     #"Fatty_Acid_zscore", 
                     #"Glutathione_zscore"
                     ), 
        scale.min = 5,
        dot.scale = 10, 
        group.by = "Age") +
  scale_color_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0
  )

fuelpath_groupedby_celltype_dotplot <- DotPlot(normscale.obj, 
                                          features = c("Glycolysis_zscore", 
                                                       "OXPHOS_zscore", 
                                                       "TCA_zscore", 
                                                       "PPP_zscore"
                                                       #, 
                                                       #"PYR_PUR_zscore", 
                                                       #"Ketone_Body_zscore", 
                                                       #"Pyruvate_zscore", 
                                                       #"FRUC_MAN_zscore", 
                                                       #"STA_SUC_zscore", 
                                                       #"PANT_COA_zscore", 
                                                       #"Fatty_Acid_zscore", 
                                                       #"Glutathione_zscore"
                                          ), 
                                          group.by = "gen_celltype") +
  scale_color_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0
  )

DoHeatmap(normscale.obj, 
        features = c("Glycolysis _zscore", 
                     "OXPHOS _zscore", 
                     "TCA _zscore"
                     #, 
                     #"PPP_zscore", 
                     #"PYR_PUR_zscore", 
                     #"Ketone_Body_zscore", 
                     #"Pyruvate_zscore", 
                     #"FRUC_MAN_zscore", 
                     #"STA_SUC_zscore", 
                     #"PANT_COA_zscore", 
                     #"Fatty_Acid_zscore", 
                     #"Glutathione_zscore"
                     ), 
        group.by = "gen_celltype")

ggsave(filename = "20260601_pretty_normscale_UMAP.pdf", plot = normscale_UMAP, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs", height = 8, width = 11, dpi = 600)

ggsave(filename = "20260601_celltype_dotplot.pdf", plot = fuelpath_groupedby_celltype_dotplot, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs", height = 8, width = 11, dpi = 600)

ggsave(filename = "20260601_age_OXPHOS_dotplot.pdf", plot = fuelpath_groupedby_age_dotplot, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs", height = 9, width = 7, dpi = 600)

ggsave(filename = "20260601_hk2_ftplot.pdf", plot = HK2_ftplot, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs", height = 8, width = 11, dpi = 600)
