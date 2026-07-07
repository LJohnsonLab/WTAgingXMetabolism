########################################################################
# Name: CCL_ST_MCC_SpaNorm_RPCA_WTAgingXMet_09.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 06. Metabolic Pathway z scoring (both global and celltype independent)  and mean normalization expression on all cells in SpaNorm RPCA obj
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
options(future.globals.maxSize = 400 * 1024^3)  # slightly under your --mem

# set seed
set.seed(42)

# to prevent this error: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so.
options(future.rng.onMisuse = "resolve")
options(future.seed = TRUE)


######### Paths ############
input_obj <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260616_SpaNorm_RPCA.qs2" 
output_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/"
date <- "20260617_" 
###########################

# NOTES:
## - prelim metabolic pathway scoring figures from 6/14 used GLOBAL z-scoring and mean expression (all cells together, NOT by cell type), done on input obj 20260612_SpaNorm_RPCA_cell_region_ann.qs2 (which did not yet have the fine_celltype col)
##       to create 20260613_SpaNorm_RPCA_scored.qs2
## - I reorganized my scripts so that I make the final SpaNorm RPCA UMAPs and imagedimplots, and calculate the cell proportions (and make those figures) before metabolically scoring 
## - I will be remaking the metabolic pathway scoring figures anyway but just wanted to put this note here 

# 01. Read in obj ----
SpaNorm_RPCA.obj <- qs_read(input_obj)

# 02. Metabolic pathway reference tables ----

## 02a. Fuel pathways reference table

met_path_table <- list(
  
  # ----fuel pathways----
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
  "Glutathione" = c("Gclc", "Gclm", "Lap3", "Gstm1", "Idh1"), 
  
  # ----amino acid pathways----
  "Ala_Asp_Glu" = c("Got1", "Gpt2", "Glud1", "Got2", "Folh1", "Aspa"), 
  "Gly_Ser_Thr" = c("Phgdh", "Dld", "Gatm", "Shmt2", "Pgam1", "Bpgm"), 
  "Cys_Met" = c("Bckdha", "Hibadh", "Auh", "Cbs", "Acad8"), 
  "Lys" = c("Ogdh", "Hadha", "Echs1"), 
  "Arg_Pro" = c("Odc1", "Sat1", "Prodh", "P4ha1", "Pycr2"),
  "His" = c("Aldh7a1", "Aldh9a1", "Cndp2"),
  "Tyr" = c("Adh1", "Adh5"),
  "Phe" = c("Maob", "Aldh3b1", "Pah"),
  "Tryp" = c("Ogdhl", "Cat", "Tdo2", "Ido1", "Tph2"),
  
  # ----lipid-related----
  "Glycerophospholipid" = c("Gpat4", "Lpgat1", "Lpin1", "Dgkd", "Pla2g15", "Pla2g7", "Ptdss1", "Lpcat2", "Etnk1", "Plpp3", "Plpp4"),
  "Sphingolipid" = c("Sptlc2", "Asah1", "Acer3", "Degs1", "Sgms1", "Kdsr", "Cers2", "Smpd1", "Cerk", "Sgpl1", "Neu1", "Sphk2"),
  
  # ----glycans----
  "Mannose_O_Glycan" = c("Fut9", "Pomgnt1", "B4galt1"),
  "N_Glycan" = c("Dpm1", "Stt3a", "Dad1", "Ddost", "Alg5", "Rpn1", "Man2a2"),
  "Mucin_O_Glycan" = c("Galnt1", "C1galt1c1"),
  "Lacto_Neolacto" = c("B3galt5", "St3gal6"),
  "Globo_Isoglobo" = c("St3gal2", "Naga"),
  "Ganglio" = c("Slc33a1", "St6galnac3", "St6galnac4"),
  "Chondroitin_Dermatan_Sulfate" = c("B3gat3", "Chsy1", "Csgalnact1", "Chst11"),
  "Heparan_Sulfate_Heparin" = c("Extl3", "Hs2st1", "Ext2"),
  "Keratan_Sulfate" = c("Fut8", "Chst1"),
  "Glycosaminoglycan_Degradation" = c("Arsb", "Gns", "Hexa", "Hexb", "Gusb", "Ids"),
  "Steroids" = c("Msmo1", "Soat1", "Comt", "Fdft1", "Ebp", "Sc5d"),
  "Arachidonic_Acid" = c("Pla2g12a", "Ptges3", "Lta4h", "Gpx1", "Ptgs1", "Cbr1", "Ltc4s"),
  "Linoleic_Acid" = c("Fads2", "Acaa1a"), 
  "D_Glutamine_Glutamate" = c("Gls")
  
)


 

library(Matrix)
# 03. GLOBAL METABOLIC PATHWAY Z-SCORING full obj  ----
pathway_zscores <- imap(met_path_table, function(genes, pathway_name) {
    
  # Step 1: get genes present in object
  pathway_genes <- intersect(genes, rownames(SpaNorm_RPCA.obj))
  
  # Step 2: extract sparse matrix (genes x cells) - no dense conversion
  mat <- GetAssayData(SpaNorm_RPCA.obj, layer = "data")[pathway_genes, , drop = F]
  
  # Step 3: z-score each gene across cells (row-wise), keeping sparse structure
  # Scale each gene: subtract mean, divide by sd
  row_means <- rowMeans(mat)
  row_sds   <- apply(mat, 1, sd)
  row_sds[row_sds == 0] <- 1  # avoid divide-by-zero for zero-variance genes
  
  
  
  colMeans((mat - row_means) / row_sds)  # broadcasts correctly
  
}) |> 
  compact() |> # drops nulls from skipped pathways 
  as.data.frame() |> 
  rename_with(~ paste0(.x, "_global_zscore")) # append zscore to pathway cols names 


# add pathway z scores for each cell to metadata 
SpaNorm_RPCA.obj <- AddMetaData(SpaNorm_RPCA.obj, metadata = pathway_zscores)

# CELL-TYPE SPECIFIC METABOLIC PATHWAY Z SCORING full obj ----

celltype_zscore_results <- map(names(met_path_table), function(pathway_name) {
  genes <- met_path_table[[pathway_name]]
  pathway_genes <- intersect(genes, rownames(SpaNorm_RPCA.obj))
  
  
  mat <- GetAssayData(SpaNorm_RPCA.obj, layer = "data")[pathway_genes, , drop = FALSE]
  
  # Get cell type for each spot
  celltypes <- SpaNorm_RPCA.obj$fine_celltype
  
  # Z-score within each cell type independently
  scores <- map(unique(celltypes), function(ct) {
    ct_cells <- names(celltypes[celltypes == ct])
    ct_mat   <- mat[, ct_cells, drop = FALSE]
    
    row_means <- rowMeans(ct_mat)
    row_sds   <- apply(ct_mat, 1, sd)
    row_sds[row_sds == 0] <- 1
    
    colMeans((ct_mat - row_means) / row_sds)
  }) %>%
    unlist()
  
  # Reorder to match original cell order
  scores[colnames(mat)]
}) %>%
  set_names(names(met_path_table)) %>%
  compact() %>%
  as.data.frame() %>%
  rename_with(~ paste0(.x, "_celltype_zscore"))

SpaNorm_RPCA.obj <- AddMetaData(SpaNorm_RPCA.obj, metadata = celltype_zscore_results)


# METABOLIC SCORING MEAN NORMALIZED EXPRESSION full obj ---- 
# Average log-normalised expression per pathway (no z-scoring)
pathway_scores_meanexp <- imap(met_path_table, function(genes, pathway_name) {
  pathway_genes <- intersect(genes, rownames(SpaNorm_RPCA.obj))
  if (length(pathway_genes) < 2) return(NULL)
  
  mat <- GetAssayData(SpaNorm_RPCA.obj, layer = "data")[pathway_genes, , drop = FALSE]
  colMeans(mat)  # average normalised expression, no scaling
}) %>%
  compact() %>%
  as.data.frame() %>%
  rename_with(~ paste0(.x, "_global_meanexp"))

SpaNorm_RPCA.obj <- AddMetaData(SpaNorm_RPCA.obj, metadata = pathway_scores_meanexp)

# save obj with metabolic scores metadata cols
qs_save(SpaNorm_RPCA.obj, paste0(output_path, date, "SpaNorm_RPCA_scored.qs2"))


# save matrices for figure making on local computer ----
## read obj back in ----
SpaNorm_RPCA.obj<- qs_read(paste0(output_path, "20260617_SpaNorm_RPCA_scored.qs2"))

## EXTRACT ALL ZSCORE METADATA AND GROUPING VARIABLES ----
zscore_cols <- grep("_zscore", colnames(SpaNorm_RPCA.obj@meta.data), value = TRUE)

# extract metadata with all zscore cols and grouping variables to load locally 
zscore_metadata <- SpaNorm_RPCA.obj@meta.data %>%
  select(gen_celltype, fine_celltype, Age, sample_ID, all_of(zscore_cols))

# save ALL zscore (global and cell specific) metadata to load locally 
qs_save(zscore_metadata, paste0(output_path, date, "SpaNorm_RPCA_zscore_metadata.qs2"))


# EXTRACT ALL MEANEXP METADATA AND GROUPING VARIABLES ----
# 1. Extract meanexp columns and age group from metadata
meanexp_cols <- grep("_meanexp", colnames(SpaNorm_RPCA.obj@meta.data), value = TRUE)


meanexp_metadata <- SpaNorm_RPCA.obj@meta.data %>%
  select(gen_celltype, fine_celltype, Age, sample_ID, all_of(meanexp_cols))

# save ALL zscore (global and cell specific) metadata to load locally 
qs_save(meanexp_metadata, paste0(output_path, date, "SpaNorm_RPCA_meanexp_metadata.qs2"))

