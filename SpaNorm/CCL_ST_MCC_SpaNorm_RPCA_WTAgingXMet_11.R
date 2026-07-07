########################################################################
# Name: CCL_ST_MCC_SpaNorm_RPCA_WTAgingXMet_11.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 11. subsetting, subclustering, and annotating astrocytes  
# Input Files: 
# Final Output Files: 
# Date created: 6/18/26
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
input_obj <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/20260617_SpaNorm_RPCA_scored.qs2" # UPDATE OBJ
output_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/"
figs_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/SpaNorm_RPCA_figs/astro_only/"
cluster_markers <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/output/SpaNorm_cluster_markers/"
genes_to_remove_path <- "/mnt/gpfs3_amd/pscratch/lajo247_uksr/WTAgingXMet/Genes_removed_from_astros.xlsx"
date <- "20260623_"
###########################



# 01. read in obj ----
SpaNorm_RPCA.obj <- qs_read(input_obj)

# 02. subset to astros ----
astro.obj <- subset(SpaNorm_RPCA.obj, subset = gen_celltype == "Astrocyte")

# 03. save astro obj to load locally ----
qs_save(astro.obj, paste0(output_path, date, "astro_SpaNorm_RPCA.qs2"))

## 03a. read in astro only obj ----
astro.obj <- qs_read(paste0(output_path, "20260618_astro_SpaNorm_RPCA.qs2"))


# SANITY CHECK
#dim(astro.obj)
#  480 194195  


UMAP_1 <- DimPlot(astro.obj)

ggsave(filename = paste0(date, "astro_notprocessed_UMAP.png"), plot = UMAP_1, path = figs_path) 

# 04. process ----
# check default assay and adjust as needed 
DefaultAssay(astro.obj) <- "Xenium"

## 04a. normalize
astro.obj <- NormalizeData(astro.obj)

## 04b. find variable features 
astro.obj <- FindVariableFeatures(astro.obj)

## 04c. scale 
astro.obj <- ScaleData(astro.obj)

## 04d. PCA
astro.obj <- RunPCA(astro.obj, npcs = 40, features = rownames(astro.obj))

PCA_age <- PCAPlot(astro.obj, group.by = "Age")

ggsave(filename = paste0(date, "astro_PCA_byage.png"), plot = PCA_age, path = figs_path) 

# elbow plot to determine number of PCs to use
elbow_plot <- ElbowPlot(astro.obj, ndims = 40, reduction = "pca")

ggsave(filename = paste0(date, "astro_elbowplot.png"), plot = elbow_plot, path = figs_path) 

# save intermediate obj ----
qs_save(astro.obj, paste0(output_path, date, "SpaNorm_RPCA_astro_nodimreduc.qs2"))

# read in intermediate obj ----
astro.obj <- qs_read(paste0(output_path, "20260619_SpaNorm_RPCA_astro_nodimreduc.qs2"))

## 03e. RunUMAP ----
astro.obj <- RunUMAP(astro.obj, dims = 1:20)  

## 03f. FindNeighbors ----
astro.obj <- FindNeighbors(astro.obj, reduction = "pca", dims = 1:20)

## 03g. FindClusters ----
astro.obj <- FindClusters(astro.obj, resolution = 0.3, graph.name = "Xenium_snn")  

# save clusters to metadata ----
Idents(astro.obj) <- "seurat_clusters"

astro.obj$clusters_res03 <- Idents(astro.obj)

## 03h. Visualize UMAP with clusters ----
UMAP_res03 <- DimPlot(astro.obj, label = T, cols = "polychrome")

imgdim_splitclusters_res03 <- ImageDimPlot(astro.obj, split.by = "seurat_clusters", fov = "fov")

imgdim_res03 <- ImageDimPlot(astro.obj, group.by = "seurat_clusters", cols = "polychrome")

# join layers before finding markers
astro.obj <- JoinLayers(astro.obj)

## 03i. find markers ----
astro.markers <- FindAllMarkers(astro.obj, group.by = "seurat_clusters", only.pos = T, logfc.threshold = 0.25)

astro.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_astro

write_csv(top5_astro, paste0(cluster_markers, date, "SpaNorm_RPCA_astro_subclusters_res03.csv"))  


# 04. feature plots for annotating ----

# astro markers 
astro_fts <- FeaturePlot(astro.obj, features = c("Aldh1l1", "Aqp4", "Gfap", "Slc7a10"))

# other cell type markers
microBAM_fts <- FeaturePlot(astro.obj, features = c("Mrc1", "P2ry12"))

oligoOPC_fts <- FeaturePlot(astro.obj, features = c("Ermn", "Pdgfra"))
neu_fts <- FeaturePlot(astro.obj, features = c("Grin2a", "Slc17a6", "Slc17a7", "Cnr1"))

# vascular
CPepend_fts <- FeaturePlot(astro.obj, features = c("Car12", "Ccdc153"))
endo_fts <- FeaturePlot(astro.obj, features = c("Flt1", "Cldn5", "Emcn"))
fibro_fts <- FeaturePlot(astro.obj, features = "Col1a2")
peri_fts <- FeaturePlot(astro.obj, features = c("Vtn", "Kcnj8", "Atp13a5"))
VLMC_fts <- FeaturePlot(astro.obj, features = c("Slc47a1", "Mgp", "Bgn"))
VSMC_fts <- FeaturePlot(astro.obj, features = c("Acta2", "Tagln"))

# save figures 
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_UMAPres03.png"), plot = UMAP_res03, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_imgdimplot_split_res03.png"), plot = imgdim_splitclusters_res03, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_imgdim_res03.png"), plot = imgdim_res03, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_ftplot.png"), plot = astro_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_microBAM_ftplot.png"), plot = microBAM_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_oligoOPC_ftplot.png"), plot = oligoOPC_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_neu_ftplot.png"), plot = neu_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_CPepend_ftplot.png"), plot = CPepend_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_endo_ftplot.png"), plot = endo_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_fibro_ftplot.png"), plot = fibro_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_peri_ftplot.png"), plot = peri_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_VLMC_ftplot.png"), plot = VLMC_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_VSMC_ftplot.png"), plot = VSMC_fts, path = figs_path)



# save fully processed obj with cluster info 
qs_save(astro.obj, paste0(output_path, date, "SpaNorm_RPCA_astro_wclus.qs2"))


# read in obj
astro.obj <- qs_read(paste0(output_path, "20260620_SpaNorm_RPCA_astro_wclus.qs2"))

#table(astro.obj$clusters_res03)
#    0     1     2     3     4     5     6     7
# 48428 43527 29585 24765 16745 12738 10605  7802


Idents(astro.obj) <- "clusters_res03"

# subset out clusters 5 and 6 bc they are contamination
astro.obj <- subset(astro.obj, idents = c("5", "6"), invert = T)

#dim(astro.obj)
#480 170852

qs_save(astro.obj, paste0(output_path, date, "SpaNorm_RPCA_astro_removed56.qs2"))
# read obj back in 
astro.obj <- qs_read(paste0(output_path, "20260621_SpaNorm_RPCA_astro_removed56.qs2"))

# remove genes
genes_to_remove <- data.frame(read_excel(genes_to_remove_path))

# create list of genes to keep
genes_to_keep <- setdiff(rownames(astro.obj), genes_to_remove$Genes)

# subset obj to genes to keep
astro.obj <- subset(astro.obj, features = genes_to_keep)

# SANITY CHECK
#dim(astro.obj)
# 377 170852

# re-process 
DefaultAssay(astro.obj) <- "Xenium" 

## Normalize ----
astro.obj <- NormalizeData(astro.obj)

## FindVariableFeatures ----
astro.obj <- FindVariableFeatures(astro.obj)

## ScaleData ----
astro.obj <- ScaleData(astro.obj)

## RunPCA ----
astro.obj <- RunPCA(astro.obj, npcs = 20, features = rownames(astro.obj))

## RunUMAP ----
## 07e. RunUMAP ----
astro.obj <- RunUMAP(astro.obj, dims = 1:20)

## 07f. FindNeighbors ----
astro.obj <- FindNeighbors(astro.obj, reduction = "pca", dims = 1:20)

## 07g. FindClusters ----
astro.obj <- FindClusters(astro.obj, resolution = 0.3, graph.name = "Xenium_snn")

# save updated res03 clusters to metadata ----
Idents(astro.obj) <- "seurat_clusters"

astro.obj$final_clusters_res03 <- Idents(astro.obj)

# find markers and save to csv file 
astro.markers_v1 <- FindAllMarkers(astro.obj, group.by = "final_clusters_res03", only.pos = T, logfc.threshold = 0.25)

astro.markers_v1 |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  ungroup() -> top10_astros_v1

write_csv(top10_astros_v1, paste0(cluster_markers, date, "top10_res03_SpaNorm_RPCA_genesremovedastro.csv"))

# dim plots 

UMAP_res03_v1 <- DimPlot(astro.obj, group.by = "final_clusters_res03", cols = "polychrome")

imgdim_clusters_res03_v1 <- ImageDimPlot(astro.obj, split.by = "final_clusters_res03", cols = "polychrome")
imgdim_res03_v1 <- ImageDimPlot(astro.obj, group.by = "final_clusters_res03", cols = "polychrome")


# ft plots 
# NOTE: astro subcluster annotations found in .../annotation_ppts/20260622_SpaNorm_RPCA_astro_annotation.pptx

astro_fts <- FeaturePlot(astro.obj, features = c("Aldh1l1", "Aqp4", "Gfap", "Slc7a10"))

LARA_fts <- FeaturePlot(astro.obj, features = c("Gja1", "Cpe", "Ubc", "Tubb2b"))
LARA_fts_2 <- FeaturePlot(astro.obj, features = c("Atp1b2", "Ptgds", "Ckb", "Psap"))
LARA_fts_3 <- FeaturePlot(astro.obj, features = c("Scg3", "Htra1", "Wif1"))

lipid_related <- FeaturePlot(astro.obj, features = c("Plin2", "Apoe", "Lpl", "Hif1a"))
lipid_related_2 <- FeaturePlot(astro.obj, features = c("App", "Clu"))

ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_cleaned_UMAP_res03.png"), plot = UMAP_res03_v1, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_cleaned_imgdimsplit_res03.png"), plot = imgdim_clusters_res03_v1, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_imgdim_res03.png"), plot = imgdim_res03_v1, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_ftplot.png"), plot = astro_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_LARA_ftplot.png"), plot = LARA_fts, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_LARA_2_ftplot.png"), plot = LARA_fts_2, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_LARA_3_ftplot.png"), plot = LARA_fts_3, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_lipid_related_ftplot.png"), plot = lipid_related, path = figs_path)
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_lipid_related_2_ftplot.png"), plot = lipid_related_2, path = figs_path)

# save final obj
qs_save(astro.obj, paste0(output_path, date, "SpaNorm_RPCA_astro_updated.qs2"))

# read object back in to get cell counts per cluster 
astro.obj <- qs_read(paste0(output_path, "20260622_SpaNorm_RPCA_astro_updated.qs2"))

#table(astro.obj$final_clusters_res03)
#0     1     2     3     4     5     6     7
#44377 30964 25997 23901 23657 13779  7781   396

# looking further at these clusters split by age for help with annotation 

#table(astro.obj$final_clusters_res03, astro.obj$Age)
#     16    36    59    92
# 0 11705 11120 11374 10178
# 1  8053  7660  7693  7558
# 2  7428  6489  5703  6377
# 3  5514  6561  6384  5442
# 4  5806  6599  5453  5799
# 5  3251  3407  3388  3733
# 6  1733  2041  1935  2072
# 7    33    45   179   139


UMAP_res03_splitage <- DimPlot(astro.obj, group.by = "final_clusters_res03", split.by = "Age",  cols = "polychrome")
ggsave(filename = paste0(date, "SpaNorm_RPCA_astro_cleaned_UMAP_res03_splitage.png"), plot = UMAP_res03_splitage, path = figs_path)

# rename idents and add to metadata ----
# make sure idents are correct
Idents(astro.obj) <- "final_clusters_res03"

astro.obj <- RenameIdents(astro.obj, 
             "0" = "Homeostatic gray matter astrocyte", 
             "1" = "Activity-responsive astrocyte", 
             "2" = "Mature homeostatic astrocyte", 
             "3" = "Metabolic homeostatic astrocyte", 
             "4" = "Cortical plate (immune-adjacent) astrocyte", 
             "5" = "Pan-reactive astrocyte", 
             "6" = "Antigen-presenting astrocyte", 
             "7" = "Interferon-responsive astrocyte"
             )


astro.obj$astro_subtype <- Idents(astro.obj)

# save object annotated obj ----
qs_save(astro.obj, paste0(output_path, date, "SpaNorm_RPCA_astro_annotated.qs2"))
