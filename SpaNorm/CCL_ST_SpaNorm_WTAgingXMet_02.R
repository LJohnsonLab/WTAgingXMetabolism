########################################################################
# Name: CCL_ST_SpaNorm_WTAgingXMet_02.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - SpaNorm 02.  running SpaNorm per sample, attaching SpaNorm assay and creating slide obj list, attaching merged SpaNorm assay to merged obj
#             - 01 
#             - 02 
#             - 03 
#             - 04 
#             - 05 
#             - 06  
#             - 07 
#             - 08  
# Input Files: 20260420_mergedobjv1_01.qs2  
# Final Output Files: 
# Date created: 6/2/26
# Last updated: 6/5/26
# Author: Chloe Lucido
########################################################################

# install SingleR
BiocManager::install("SingleR")  
BiocManager::install("bhuvad/SpaNorm")


# Load Libraries ----
library(spacexr)
library(Seurat)
library(SeuratDisk)
library(future)
library(ggplot2)
library(arrow)
library(hdf5r) # for anndata -> seurat conversion
library(presto)
library(glmGamPoi)
library(readr)
library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)
library(patchwork)
library(sceasy) # for anndata -> seurat conversion
library(reticulate) # for anndata -> seurat conversion
library(SPLIT)
library(qs2) 
library(RColorBrewer)
library(Polychrome)
library(purrr)
library(readxl)
library(BPCells)
library(SpaNorm)
library(SpatialExperiment)

options(future.globals.maxSize = 100 *1024^3)

# set seed
set.seed(42)

####### Paths and Obj names #######
inputobj_path <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260420_mergedobjv1_01.qs2"
output_path <- "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/"
output_WTmergedobj_name <- "20260604_mergedobj_SpaNorm_Assay_noscaling.qs2"
output_SpaNorm_matrices_list_name <- "20260604_samplelist_of_SpaNorm_matrices.qs2"
  
# 01. read in merged obj ----
WTmerged.obj <- qs_read(inputobj_path)


## 01a. Update merged obj ----
WTmerged.obj <- UpdateSeuratObject(WTmerged.obj)

# 02. map through merged obj and subset by sample (Adapted from Jose Arbones-Mainar)----


#######################################################################
# FUNCTION: subset mergedobj into one sample (using sample ID), attach coordinates from corresponding fov (utilizing fov_lookup), convert seurat obj to spatial experiment obj to then run SpaNorm
spanorm_one_sample <- function(s, xen, fov_lookup) {
  # s = sample ID; xen = xenium seurat obj; fov_lookup = table corresponding fov name to slide ID
  
  # subset to only one sample from merged obj
  seu_s <- subset(xen, subset = sample_ID == s)
  
  # get slide ID information from sample metadata
  sample_meta <- xen@meta.data |>
    filter(sample_ID == s) |>
    slice(1) |>
    select(slide)
  
  slide_s <- sample_meta$slide 
  
  # extract coords for SpatialExperiment (spe) object
  # fov_lookup maps sample ID to FOV name , after subsetting the fov image obj is retained but the list of cells needs to be updated
  coords <- GetTissueCoordinates(seu_s, image = fov_lookup[[s]]) |>
    filter(cell %in% Cells(seu_s))
  
  # build spatial experiment object for each sample in order to run SpaNorm
  spe <- SpatialExperiment(
    assays = list(counts = GetAssayData(seu_s, assay = "Xenium", layer = "counts")),
    spatialCoords = as.matrix(coords[, c("x", "y")])
  )
  
  # now run SpaNorm on each individual sample obj
  spe <- SpaNorm(spe,
                 sample.p = 0.05,
                 df.tps = 6,
                 adj.method = "auto",
                 backend = "auto",
                 verbose = T
  )
  # extract normalized matrix from SpaNorm and restore Seurat barcodes
  out <- SummarizedExperiment::assay(spe, "logcounts")
  
  colnames(out) <- Cells(seu_s) # resetting colnames to cell barcodes
  
  # attach sample-level metadata back to extracted matrices via attributes
  attr(out, "sample_ID") <- s
  attr(out, "slide") <- slide_s
  out
}

#########################################################################


# extract metadata from obj to be referenced later ----
samples <- WTmerged.obj@meta.data |>
  distinct(sample_ID, slide)
# Build the sample_id -> FOV-name lookup from the object itself; LoadXenium uses default name "fov" and merge renames collisions
# (so the actual names are usually fov.1 / fov.2, not the slide labels).
fov_to_slide <- map_chr(Images(WTmerged.obj), \(img) {
  # Each FOV's cells share the same slide value in meta.data; use the first one as the slide tag.
  unique(WTmerged.obj@meta.data[Cells(WTmerged.obj[[img]]), "slide"])[1]
}) |> set_names(Images(WTmerged.obj))

# Invert: slide -> FOV, then index by sample slide to get the per-sample FOV lookup.
slide_to_fov <- set_names(names(fov_to_slide), unname(fov_to_slide))

fov_lookup <- set_names(slide_to_fov[samples$slide], samples$sample_id)


# Map over sample IDs — no tribble needed anywhere.
norm_matrices <- map(
  set_names(samples$sample_ID),
  spanorm_one_sample,
  xen        = WTmerged.obj,
  fov_lookup = fov_lookup
)

# Bind the per-sample matrices into one wide normalized matrix and reorder to the merged object's cell order.
# Qualify reduce() because SpatialExperiment loads IRanges, whose reduce() masks purrr::reduce.
norm_mat <- purrr::reduce(norm_matrices, cbind)

norm_mat <- norm_mat[, Cells(WTmerged.obj)]
WTmerged.obj[["SpaNorm"]] <- CreateAssay5Object(data = norm_mat)


DefaultAssay(WTmerged.obj) <- "SpaNorm"

# save files 
qs_save(WTmerged.obj, paste0(output_path, output_WTmergedobj_name))

qs_save(SpaNorm_obj_list, paste0(output_path, output_SpaNorm_slide_objlist_name))
