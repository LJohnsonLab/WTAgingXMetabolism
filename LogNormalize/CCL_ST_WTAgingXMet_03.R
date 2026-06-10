########################################################################
# Name: CCL_ST_WTAgingXMet_03.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - 03. annotating all cells from part 2 normalized and scaled merged obj   
#             - 01 read in merged obj
#             - 02 annotate clusters
#             - 03 assigning cell types to clusters (renaming idents and storing in metadata)
#             - 04 saving final annotated obj => 20260428_mergedobj_gen_annotations_03.qs2
# Input Files: 20260424_mergedobj_res04_02.qs2  
# Final Output Files: 20260428_mergedobj_gen_annotations_03.qs2
# Date created: 4/24/26
# Last updated: 4/28/26
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

# Notes:
# - Old res 0.3 UMAP SAVED AS 20260420_allcells_UMAP_normalize_scale_res03.png in /figs/UMAPs


# 01. read in merged obj ----

# left off part 2 using 30 dims for umap and neighbors and clustering at res = 0.4

WTmerged.obj <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260424_mergedobj_res04_02.qs2")


# 02. Annotating ----

## 02a. save res 0.4 clusters (from part 2) ----

WTmerged.obj$res04_clusters <- Idents(WTmerged.obj)

## 02b. UMAP
UMAP_res04 <- DimPlot(WTmerged.obj, label = T, cols = "polychrome")

## 02c. save UMAP to reference 
ggsave(filename = "20260424_normalize_scale_UMAP_res04.png", plot = UMAP_res04, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs/", width = 9, height = 10)


## 02b. FindAllMarkers ----
WTmerged.obj.markers <- FindAllMarkers(WTmerged.obj, group.by = "res04_clusters")

WTmerged.obj.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_WTmerged

write_csv(top5_WTmerged, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/cluster_markers/20260424_top5_WT_xenium_res04.csv")

WTmerged.obj.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  ungroup() -> top10_WTmerged

write_csv(top10_WTmerged, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/cluster_markers/20260427_top10_WT_xenium_res04.csv")


# annotate clusters 

FeaturePlot(WTmerged.obj, features = c("Aldh1l1", "Gfap", "Slc7a10", "Aqp4")) # astro

FeaturePlot(WTmerged.obj, features = c("Tmem119", "Aif1", "Trem2", "P2ry12")) # mircoglia

FeaturePlot(WTmerged.obj, features = c("Cst7", "Clec7a", "Siglech"))

FeaturePlot(WTmerged.obj, features = c("Ermn", "Cldn11", "Mog")) # oligo

FeaturePlot(WTmerged.obj, features = c("Ccdc153", "Dnah11", "Tmem212")) # epend

FeaturePlot(WTmerged.obj, features = c("Pdgfra", "Tnr")) # OPC

FeaturePlot(WTmerged.obj, features = c("Mgp", "Bgn", "Slc47a1")) # VLMC

FeaturePlot(WTmerged.obj, features = c("Acta2", "Tagln")) # VSMC

FeaturePlot(WTmerged.obj, features = c("Vtn", "Kcnj8", "Atp13a5")) # pericytes

FeaturePlot(WTmerged.obj, features = c("Flt1", "Emcn", "Cldn5")) # endo

FeaturePlot(WTmerged.obj, features = c("Car12", "Ttr")) # CP

FeaturePlot(WTmerged.obj, features = c("Col1a2", "Lum")) # fibroblast

FeaturePlot(WTmerged.obj, features = c("Dcx", "Pax6")) # NPC

FeaturePlot(WTmerged.obj, features = c("Mgl2", "Mrc1")) # BAM

FeaturePlot(WTmerged.obj, features = c("Cd3d", "Gzmb", "Pdcd1")) # T cell

FeaturePlot(WTmerged.obj, features = c("Snap25", "Grin2a")) # pan-neu

FeaturePlot(WTmerged.obj, features = c("Fezf2", "Slc17a7", "Slc17a6", "Htr2c")) # general excitatory neu

FeaturePlot(WTmerged.obj, features = c("Cux2", "Otof", "Stard8", "Lypd1", "Lrg1")) # L2/3 IT (exc) neu

FeaturePlot(WTmerged.obj, features = "Rapgef3") # L5 PT (exc) neu

FeaturePlot(WTmerged.obj, features = c("Rorb", "Tcap", "Rspo1", "Whrn")) # L5 IT (exc) neu

FeaturePlot(WTmerged.obj, features = "Tunar") # L6 IT (exc) neu

FeaturePlot(WTmerged.obj, features = "Pou3f1") # L5 ET (exc) neu

FeaturePlot(WTmerged.obj, features = c("Tshz2", "Tox2")) # L5/6 NP (exc) neu

FeaturePlot(WTmerged.obj, features = c("Syt6", "Trh")) # L6 CT (exc) neu

FeaturePlot(WTmerged.obj, features = c("Gad1", "Lhx6", "Cnr1", "Serpinf1", "Slc32a1")) # general inhibitory neu

FeaturePlot(WTmerged.obj, features = c("Lamp5", "Ndnf")) # LAMP5 (inh) neu

FeaturePlot(WTmerged.obj, features = c("Pvalb", "Nts", "Tac1")) # PVALB (inh) neu

FeaturePlot(WTmerged.obj, features = "Sncg") # SNCG (inh) neu

FeaturePlot(WTmerged.obj, features = c("Sst", "Calb1")) # SST (inh) neu

FeaturePlot(WTmerged.obj, features = c("Vip", "Calb2", "Pthlh", "Crh")) # VIP (inh) neu

FeaturePlot(WTmerged.obj, features = c("Ccl5", "Nkg7", "Klrb1c", "Il1b")) # some top markers from cluster 21


imgdim_fov1_clusters <- ImageDimPlot(WTmerged.obj, split.by = "res04_clusters", cols = "polychrome")

imgdim_fov1 <- ImageDimPlot(WTmerged.obj, group.by = "res04_clusters", cols = "polychrome")

ggsave(filename = "20260427_normalize_scale_imgdim_clusters_res04.png", plot = imgdim_fov1_clusters, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/ImageDimPlot/", width = 11, height = 10)
ggsave(filename = "20260427_normalize_scale_imgdim_res04.png", plot = imgdim_fov1, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/ImageDimPlot/", width = 11, height = 10)


## 02c. finding subclusters for cluster 2 ----
# cluster 2 had a mix of GABAergic and Glutamatergic markers and regional spread so I wanted to see if I could parse these out 
# refer to 20260424_cellannotation.pptx


WTmerged.obj <- FindSubCluster(WTmerged.obj, cluster = "2", resolution = 0.2, graph.name = "Xenium_snn")

UMAP_res04_sub2 <- DimPlot(WTmerged.obj, group.by = "sub.cluster", cols = "polychrome")

imgdim_fov1_sub2clusters <- ImageDimPlot(WTmerged.obj, split.by = "sub.cluster", cols = "polychrome")

Idents(WTmerged.obj) <- "sub.cluster"

# find markers for subclustered 2

WTmerged.obj.sub2markers <- FindAllMarkers(WTmerged.obj, group.by = Idents(WTmerged.obj))

WTmerged.obj.sub2markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_WTmerged_sub2

write_csv(top5_WTmerged_sub2, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/cluster_markers/20260427_top5_WTxenium_res04_sub2.csv")

# save subclusters to metadata 

WTmerged.obj$sub2_clusters <- Idents(WTmerged.obj)

# conclusion: still couldn't parse out GABA- vs glutamatergic neurons, will subset all the neurons together and hopefully I will get a clear picture

## 02d. finding subclusters for cluster 15 ----
# cluster 15 contained a mix of oligodendrocyte markers and DAM markers, when feature plotted for Oligo markers, the bottom half of the cluster highly expressed the markers but the other half did not 
# refer to 20260424_cellannotation.pptx

WTmerged.obj <- FindSubCluster(WTmerged.obj, cluster = "15", resolution = 0.1, graph.name = "Xenium_snn")

UMAP_res04_sub15 <- DimPlot(WTmerged.obj, group.by = "sub.cluster", cols = "polychrome")

imgdim_fov1_sub15clusters <- ImageDimPlot(WTmerged.obj, split.by = "sub.cluster", cols = "polychrome")

Idents(WTmerged.obj) <- "sub.cluster"

WTmerged.obj.sub15markers <- FindAllMarkers(WTmerged.obj, group.by = Idents(WTmerged.obj))

WTmerged.obj.sub15markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_WTmerged_sub15

write_csv(top5_WTmerged_sub15, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/cluster_markers/20260427_top5_WTxenium_res04_sub15.csv")

# save subclusters in metadata 
WTmerged.obj$sub2_15_clusters <- Idents(WTmerged.obj)


# save UMAP and imagedimplot with both 2 and 15 subclusters 

ggsave(filename = "20260428_normalize_scale_UMAP_sub2_15_res04.png", plot = UMAP_res04_sub15, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs/", width = 11, height = 10)
ggsave(filename = "20260428_normalize_scale_imgdim_2_15_subclust_res04.png", plot = imgdim_fov1_sub15clusters, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/ImageDimPlot/", width = 11, height = 10)



# 03. Renaming clusters to reflect their cell types ----
# this is a preliminary annotation and will be the most general cell typing

# make sure the idents are set correctly 
Idents(WTmerged.obj)

# subset out cluster 22 bc it only has 2 cells 
WTmerged.obj <- subset(WTmerged.obj, idents = "22", invert = T)

# RenameIdents 

WTmerged.obj <- RenameIdents(WTmerged.obj, 
                             "0" = "Oligodendrocyte",
                             "1" = "Astrocyte", 
                             "2_0" = "Neuron", 
                             "2_1" = "Neuron", 
                             "2_2" = "Neuron", 
                             "2_3" = "Neuron", 
                             "2_4" = "Neuron",
                             "3" = "Neuron", 
                             "4" = "Neuron", 
                             "5" = "Endothelial", 
                             "6" = "Vascular Mural", 
                             "7" = "Neuron", 
                             "8" = "Microglia", 
                             "9" = "Neuron", 
                             "10" = "Neuron", 
                             "11" = "Neuron", 
                             "12" = "OPC", 
                             "13" = "Oligodendrocyte", 
                             "14" = "Neuron", 
                             "15_0" = "Oligodendrocyte", 
                             "15_1" = "Oligodendrocyte", 
                             "15_2" = "Oligodendrocyte", 
                             "16" = "Pericyte", 
                             "17" = "Ependymal", 
                             "18" = "CP", 
                             "19" = "Neuron", 
                             "20" = "Neuron", 
                             "21"= "Immune"
)

# save cell typed clusters into metadata
WTmerged.obj$gen_celltype <- Idents(WTmerged.obj)


# SANITY CHECK
gen_celltype_UMAP <- DimPlot(WTmerged.obj, cols = "polychrome")

gen_celltype_imgdimplot <- ImageDimPlot(WTmerged.obj, group.by = "gen_celltype", fov = "fov", cols = "polychrome")

gen_celltype_clusters_imgdimplot <- ImageDimPlot(WTmerged.obj, split.by = "gen_celltype", fov = "fov", cols = "polychrome")


ggsave(filename = "20260428_normalize_scale_gen_celltype_UMAP.png", plot = gen_celltype_UMAP, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/UMAPs/", width = 11, height = 10)
ggsave(filename = "20260428_normalize_scale_gen_celltype_imgdim.png", plot = gen_celltype_imgdimplot, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/ImageDimPlot/", width = 11, height = 10)
ggsave(filename = "20260428_normalize_scale_gen_celltype_clusters_imgdim.png", plot = gen_celltype_clusters_imgdimplot, path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/figs/ImageDimPlot/", width = 11, height = 10)


# 04. Save final annotated obj ----
qs_save(WTmerged.obj, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260428_mergedobj_gen_annotations_03.qs2")
