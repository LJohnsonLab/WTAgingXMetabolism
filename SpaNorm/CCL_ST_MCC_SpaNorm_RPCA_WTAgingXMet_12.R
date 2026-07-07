########################################################################
# Name: CCL_ST_MCC_SpaNorm_RPCA_WTAgingXMet_12.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 12. metabolic scoring on SpaNorm RPCA astrocyte object, making pretty figures from astrocyte obj for ICBEM 2026 poster and pres  
# Input Files: 
# Final Output Files: 
# Date created: 6/23/26
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
library(SCP) # for pretty UMAPs


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
input_obj <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260623_SpaNorm_RPCA_astro_annotated.qs2" 
output_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/"
figs_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/SpaNorm_RPCA_figs/astro_only/"
date <- "20260624_"
###########################



# 01. read in obj ----
astro.obj <- qs_read(input_obj)

# 02. reorder idents for plotting ----
## 02a. ensure correct idents ---
Idents(astro.obj) <- "astro_subtype"

## 02b. reorder idents ----
astro.obj$astro_subtype <- factor(astro.obj$astro_subtype, 
                                  levels = c(
                                           "Homeostatic gray matter astrocyte",
                                           "Metabolic homeostatic astrocyte",
                                           "Mature homeostatic astrocyte",
                                           "Cortical plate (immune-adjacent) astrocyte", 
                                           "Activity-responsive astrocyte", 
                                           "Antigen-presenting astrocyte",
                                           "Interferon-responsive astrocyte", 
                                           "Pan-reactive astrocyte"))


# 03. astro figs  ----

## 03a. UMAP ----
astro_UMAP <- CellDimPlot(
  srt = astro.obj,
  group.by = "astro_subtype",
  reduction = "umap",
  theme_use = "theme_blank",
  pt.size = 0.6) + scale_color_manual(values = RColorBrewer::brewer.pal(8, "Paired")) +
                        guides(color = guide_legend(title = NULL, # remove legend title
                                                    ncol = 1, # keep legend in 1 column
                                                    override.aes = list(size = 3))) # legend dot size

## 03b. imgdimplots ----
### grouped by cell type, not split ----
astro_imagedim_fov7_24 <- ImageDimPlot(astro.obj,
                              group.by = "astro_subtype", 
                              cols = RColorBrewer::brewer.pal(8, "Paired"),
                              dark.background = F,
                              size = 0.9, 
                              fov = "fov.7.24") + labs(fill = NULL)

### grouped and split ----
astro_imagedim_split_fov7_24 <- ImageDimPlot(astro.obj,
                              group.by = "astro_subtype", 
                              split.by = "astro_subtype",
                              cols = RColorBrewer::brewer.pal(8, "Paired"),
                              dark.background = F,
                              size = 0.9,
                              fov = "fov.7.24") + labs(fill = NULL)



# metabolic pathway figures ----
astro_pyruvate_zscore <- FeaturePlot(astro.obj, 
            features = "Pyruvate_global_zscore") + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

astro_lacto_neolacto_zscore <- FeaturePlot(astro.obj, 
            features = "Lacto_Neolacto_global_zscore") + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

astro_glutathione_zscore <- FeaturePlot(astro.obj, 
            features = "Glutathione_global_zscore") + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

astro_glycero_zscore <- FeaturePlot(astro.obj, 
            features = "Glycerophospholipid_global_zscore") + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

astro_sphingo_zscore <- FeaturePlot(astro.obj, 
            features = "Sphingolipid_global_zscore") + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

astro_mucin_o_glycan_zscore <- FeaturePlot(astro.obj, 
            features = "Mucin_O_Glycan_global_zscore") + ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

## 04. save these plots ----
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_subtype_UMAP.png"), plot = astro_UMAP, path = figs_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_subtype_imgdimplot.png"), plot = astro_imagedim_fov7_24, path = figs_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_subtype_imgdimplot_split.png"), plot = astro_imagedim_split_fov7_24, path = figs_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_pyruvate_ftplot.png"), plot = astro_pyruvate_zscore, path = figs_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_lacto_neolacto_ftplot.png"), plot = astro_lacto_neolacto_zscore, path = figs_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_glutathione_ftplot.png"), plot = astro_glutathione_zscore, path = figs_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNorm_RPCA_glycerophospholipid_ftplot.png"), plot = astro_glycero_zscore, path = figs_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_sphingolipid_ftplot.png"), plot = astro_sphingo_zscore, path = figs_path, width = 11, height = 10)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_mucin_o_glycan_ftplot.png"), plot = astro_mucin_o_glycan_zscore, path = figs_path, width = 11, height = 10)


# 05. Metabolic pathway reference tables ----

## 05a. Fuel pathways reference table

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


 

#library(Matrix)
# 03. GLOBAL METABOLIC PATHWAY Z-SCORING full obj  ----
pathway_zscores <- imap(met_path_table, function(genes, pathway_name) {
    
  # Step 1: get genes present in object
  pathway_genes <- intersect(genes, rownames(astro.obj))
  
  # Step 2: extract sparse matrix (genes x cells) - no dense conversion
  mat <- GetAssayData(astro.obj, layer = "data")[pathway_genes, , drop = F]
  
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
astro.obj <- AddMetaData(astro.obj, metadata = pathway_zscores)

# CELL-TYPE SPECIFIC METABOLIC PATHWAY Z SCORING full obj ----

celltype_zscore_results <- map(names(met_path_table), function(pathway_name) {
  genes <- met_path_table[[pathway_name]]
  pathway_genes <- intersect(genes, rownames(astro.obj))
  
  
  mat <- GetAssayData(astro.obj, layer = "data")[pathway_genes, , drop = FALSE]
  
  # Get cell type for each spot
  celltypes <- astro.obj$astro_subtype
  
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
  rename_with(~ paste0(.x, "_subtype_zscore"))


# add subtype independent pathway z scores
astro.obj <- AddMetaData(astro.obj, metadata = celltype_zscore_results)

# save obj with metabolic scores metadata cols
qs_save(astro.obj, paste0(output_path, date, "SpaNorm_RPCA_astro_scored.qs2"))


# save matrices for figure making on local computer ----

## EXTRACT ALL ZSCORE METADATA AND GROUPING VARIABLES ----
zscore_cols <- grep("_zscore", colnames(astro.obj@meta.data), value = TRUE)

# extract metadata with all zscore cols and grouping variables to load locally 
zscore_metadata <- astro.obj@meta.data %>%
  select(astro_subtype, Age, sample_ID, all_of(zscore_cols))

# save ALL zscore (global and cell specific) metadata to load locally 
qs_save(zscore_metadata, paste0(output_path, date, "SpaNorm_RPCA_astro_zscore_metadata.qs2"))

