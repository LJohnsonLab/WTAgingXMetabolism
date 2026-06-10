########################################################################
# Name: CCL_ST_SpaNorm_RPCA_WTAgingXMet_04.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - SpaNorm 04.  annotating SpaNorm RPCA integrated obj
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
# Date created: 6/8/26
# Last updated: 6/8/26
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


#########Paths##########
root_path_to_obj <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260605_SpaNorm_RPCAint_processed.qs2"
fig_output_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/SpaNorm_RPCA_figs/"


# read in obj ----
SpaNorm_RPCA.obj <- qs_read(root_path_to_obj)


UMAP <- DimPlot(SpaNorm_RPCA.obj, group.by = "clusters_res05", label = T, cols = "polychrome")

ImgUMAP <- ImageDimPlot(SpaNorm_RPCA.obj, split.by = "clusters_res05", cols = "polychrome")


# annotate clusters 

astro_ftplot <- FeaturePlot(WTmerged.obj, features = c("Aldh1l1", "Gfap", "Slc7a10", "Aqp4")) # astro

microglia_ftplot <- FeaturePlot(WTmerged.obj, features = c("Tmem119", "Aif1", "Trem2", "P2ry12")) # mircoglia

microglia1_ftplot <- FeaturePlot(WTmerged.obj, features = c("Cst7", "Clec7a", "Siglech"))

oligo_ftplot <- FeaturePlot(WTmerged.obj, features = c("Ermn", "Cldn11", "Mog")) # oligo

epend_ftplot <- FeaturePlot(WTmerged.obj, features = c("Ccdc153", "Dnah11", "Tmem212")) # epend

OPC_ftplot <- FeaturePlot(WTmerged.obj, features = c("Pdgfra", "Tnr")) # OPC

VLMC_ftplot <- FeaturePlot(WTmerged.obj, features = c("Mgp", "Bgn", "Slc47a1")) # VLMC

VSMC_ftplot <- FeaturePlot(WTmerged.obj, features = c("Acta2", "Tagln")) # VSMC

peri_ftplot <- FeaturePlot(WTmerged.obj, features = c("Vtn", "Kcnj8", "Atp13a5")) # pericytes

endo_ftplot <- FeaturePlot(WTmerged.obj, features = c("Flt1", "Emcn", "Cldn5")) # endo

CP_ftplot <- FeaturePlot(WTmerged.obj, features = c("Car12", "Ttr")) # CP

fibroblast_ftplot <- FeaturePlot(WTmerged.obj, features = c("Col1a2", "Lum")) # fibroblast

NPC_ftplot <- FeaturePlot(WTmerged.obj, features = c("Dcx", "Pax6")) # NPC

BAM_ftplot <- FeaturePlot(WTmerged.obj, features = c("Mgl2", "Mrc1")) # BAM

Tcell_ftplot <- FeaturePlot(WTmerged.obj, features = c("Cd3d", "Gzmb", "Pdcd1")) # T cell

panneu_ftplot <- FeaturePlot(WTmerged.obj, features = c("Snap25", "Grin2a")) # pan-neu

excneu_ftplot <- FeaturePlot(WTmerged.obj, features = c("Fezf2", "Slc17a7", "Slc17a6", "Htr2c")) # general excitatory neu

L23_IT_ftplot <- FeaturePlot(WTmerged.obj, features = c("Cux2", "Otof", "Stard8", "Lypd1", "Lrg1")) # L2/3 IT (exc) neu

L5_PT_ftplot <- FeaturePlot(WTmerged.obj, features = "Rapgef3") # L5 PT (exc) neu

L5_IT_ftplot <- FeaturePlot(WTmerged.obj, features = c("Rorb", "Tcap", "Rspo1", "Whrn")) # L5 IT (exc) neu

L6_IT_ftplot <- FeaturePlot(WTmerged.obj, features = "Tunar") # L6 IT (exc) neu

L5_ET_ftplot <- FeaturePlot(WTmerged.obj, features = "Pou3f1") # L5 ET (exc) neu

L56_NP_ftplot <- FeaturePlot(WTmerged.obj, features = c("Tshz2", "Tox2")) # L5/6 NP (exc) neu

L6_CT_ftplot <- FeaturePlot(WTmerged.obj, features = c("Syt6", "Trh")) # L6 CT (exc) neu

inhneu_ftplot <- FeaturePlot(WTmerged.obj, features = c("Gad1", "Lhx6", "Cnr1", "Serpinf1", "Slc32a1")) # general inhibitory neu

lamp5_ftplot <- FeaturePlot(WTmerged.obj, features = c("Lamp5", "Ndnf")) # LAMP5 (inh) neu

pvalb_ftplot <- FeaturePlot(WTmerged.obj, features = c("Pvalb", "Nts", "Tac1")) # PVALB (inh) neu

sncg_ftplot <- FeaturePlot(WTmerged.obj, features = "Sncg") # SNCG (inh) neu

sst_ftplot <- FeaturePlot(WTmerged.obj, features = c("Sst", "Calb1")) # SST (inh) neu

vip_ftplot <- FeaturePlot(WTmerged.obj, features = c("Vip", "Calb2", "Pthlh", "Crh")) # VIP (inh) neu


# save plots 

ggsave(filename = "20260608_SpaNormRPCA_UMAP.png", plot = UMAP, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_imgdimplot.png", plot = ImgUMAP, path = fig_output_path, width = 11, height = 10)

ggsave(filename = "20260608_SpaNormRPCA_astro_ftplot.png", plot = astro_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_micro_ftplot.png", plot = microglia_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_micro1_ftplot.png", plot = microglia1_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_oligo_ftplot.png", plot = oligo_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_epend_ftplot.png", plot = epend_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_OPC_ftplot.png", plot = OPC_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_VLMC_ftplot.png", plot = VLMC_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_VSMC_ftplot.png", plot = VSMC_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_peri_ftplot.png", plot = peri_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_endo_ftplot.png", plot = endo_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_CP_ftplot.png", plot = CP_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_fibroblast_ftplot.png", plot = fibroblast_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_NPC_ftplot.png", plot = NPC_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_BAM_ftplot.png", plot = BAM_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_Tcell_ftplot.png", plot = Tcell_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_panneu_ftplot.png", plot = panneu_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_excneu_ftplot.png", plot = excneu_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_L23IT_ftplot.png", plot = L23_IT_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_L5PT_ftplot.png", plot = L5_PT_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_L5IT_ftplot.png", plot = L5_IT_ftplot, path = fig_output_path, width = 11, height = 10)

ggsave(filename = "20260608_SpaNormRPCA_L6IT_ftplot.png", plot = L6_IT_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_L5ET_ftplot.png", plot = L5_ET_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_L56NP_ftplot.png", plot = L56_NP_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_L6CT_ftplot.png", plot = L6_CT_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_inhneu_ftplot.png", plot = inhneu_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_lamp5_ftplot.png", plot = lamp5_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_pvalb_ftplot.png", plot = pvalb_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_sncg_ftplot.png", plot = sncg_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_sst_ftplot.png", plot = sst_ftplot, path = fig_output_path, width = 11, height = 10)
ggsave(filename = "20260608_SpaNormRPCA_vip_ftplot.png", plot = vip_ftplot, path = fig_output_path, width = 11, height = 10)
