########################################################################
# Name: CCL_ST_WTAgingXMet_OLD_01.R
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         - v1 01. creating xenium objs, subsetting (nCount > 5), merging files, attaching metadata
#           - contains BPCells raw count matrices
# Input Files:  
# Output Files: 
# Date created: 4/15/26
# Last updated: 4/22/26
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


# 01. Loading in Xenium data ----

### 01a. Run 1 Paths ----

#### slide 0021999
path1 <- "/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run1_20250516__201914__20250516_AgingXMetabolism_1/Slide1_output-XETG00118__0021999__Region_1__20250516__201923"

#### slide 0022573
path2 <- "/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run1_20250516__201914__20250516_AgingXMetabolism_1/Slide2_output-XETG00118__0022573__Region_1__20250516__201924"


### 01b. Run 2 Paths ----

#### Slide 0021991
path3 <- "/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run2_20250606__202944__20250606_AgingXMetabolism_2/Slide1_output-XETG00118__0021991__Region_1__20250606__202953"

#### Slide 0021998
path4 <- "/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run2_20250606__202944__20250606_AgingXMetabolism_2/Slide2_output-XETG00118__0021998__Region_1__20250606__202953"


### 01c. Run 4 Paths ----

#### Slide 0069039
path5 <- "/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run4_20250919__205218__20250919_AgingXMetabolism_4/Slide1_output-XETG00118__0069039__Region_1__20250919__205228"

#### Slide 0069045
path6 <- "/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run4_20250919__205218__20250919_AgingXMetabolism_4/Slide2_output-XETG00118__0069045__Region_1__20250919__205229"


### 01d. Run 5 Paths ----

#### Slide 0069118
path7 <- "/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run5_20260306_210910_20260306_AgingXMetabolism_5/slide1_output-XETG00118__0069118__Region_1__20260306__210910"


### 01e. Create slide_ID and path dictionary and create path df ----
# Named vector — names are Slide_IDs matching your metadata sheet
paths <- c(
  "0021999" = path1,
  "0022573" = path2,
  "0021991" = path3,
  "0021998" = path4,
  "0069039" = path5,
  "0069045" = path6,
  "0069118" = path7
)

# make dataframe from paths list 
paths <- c(path1, path2, path3, path4, path5, path6, path7)


# 02. Create Xenium Objs ----

## 02a. map through each path in paths df and creates xenium objs (stored in xenium_list) ----
xenium_list <- map(paths, function(p) {LoadXenium(p, fov = "fov", segmentations = "cell")})

## 02b. Rename cells using Slide_ID from metadata ----

# adding slide number to cell names (using last 3 digits of slide ID as cell ID prefixes)

xenium.obj1 <- RenameCells(xenium_list[[1]], add.cell.id = "slide999")

xenium.obj2 <- RenameCells(xenium_list[[2]], add.cell.id = "slide573")

xenium.obj3 <- RenameCells(xenium_list[[3]], add.cell.id = "slide991")

xenium.obj4 <- RenameCells(xenium_list[[4]], add.cell.id = "slide998")

xenium.obj5 <- RenameCells(xenium_list[[5]], add.cell.id = "slide039")

xenium.obj6 <- RenameCells(xenium_list[[6]], add.cell.id = "slide045")

xenium.obj7 <- RenameCells(xenium_list[[7]], add.cell.id = "slide118")



## 02c. Create metadata column to store slide IDs ----
xenium.obj1$slide <- names(paths)[1] 
xenium.obj2$slide <- names(paths)[2]
xenium.obj3$slide <- names(paths)[3]
xenium.obj4$slide <- names(paths)[4]
xenium.obj5$slide <- names(paths)[5]
xenium.obj6$slide <- names(paths)[6]
xenium.obj7$slide <- names(paths)[7]


# 03. Attach full brain XE annotations ----


## 03a. paths to each samples XE annotation ----
files <- list(
  # run 1
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run1_slide0021999_92wk_M1_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run1_slide0022573_59wk_F2_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run1_slide0022573_59wk_M1_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run1_slide0022573_92wk_F2_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run1_slide0022573_92wk_M2_cells_stats.csv",
  
  # run 2
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run2_slide0021991_16wk_M1_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run2_slide0021991_36wk_F1_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run2_slide0021991_36wk_M2_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run2_slide0021998_16wk_F1_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run2_slide0021998_16wk_M2_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run2_slide0021998_36wk_M1_cells_stats.csv",
  
  # run 4
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run4_slide0069039_16wk_M4_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run4_slide0069039_59wk_F5_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run4_slide0069039_92wk_F3_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run4_slide0069039_92wk_M3_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run4_slide0069045_16wk_F4_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run4_slide0069045_36wk_F3_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run4_slide0069045_59wk_M3_cells_stats.csv",
  
  
  # run 5
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run5_slide0069118_16wk_F2_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run5_slide0069118_36wk_F4_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run5_slide0069118_36wk_M4_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run5_slide0069118_59wk_F3_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run5_slide0069118_59wk_M5_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run5_slide0069118_92wk_F5_cells_stats.csv"
  
)

## 03b. formatting XE full brain annotations to add to seurat metadata ----

# reading and tagging each csv cell stats file with corresponding sample name 

annotation_list <- map(files, function(x) { 
  
  # reading each csv file, skipping first 2 metadata lines
  
  read.csv(x, skip = 2, blank.lines.skip = TRUE)
  
})


# defining which samples are in each xenium object for proper prefix placement

prefix_map <- c(
  "92wk_M1" = "slide999_", 
  "59wk_F2" = "slide573_", 
  "59wk_M1" = "slide573_", 
  "92wk_F2" = "slide573_", 
  "92wk_M2" = "slide573_", 
  "16wk_M1" = "slide991_", 
  "36wk_F1" = "slide991_", 
  "36wk_M2" = "slide991_", 
  "16wk_F1" = "slide998_", 
  "16wk_M2" = "slide998_", 
  "36wk_M1" = "slide998_", 
  "16wk_M4" = "slide039_", 
  "59wk_F5" = "slide039_", 
  "92wk_F3" = "slide039_", 
  "92wk_M3" = "slide039_", 
  "16wk_F4" = "slide045_", 
  "36wk_F3" = "slide045_", 
  "59wk_M3" = "slide045_", 
  "16wk_F2" = "slide118_", 
  "36wk_F4" = "slide118_",
  "36wk_M4" = "slide118_",
  "59wk_F3" = "slide118_",
  "59wk_M5" = "slide118_",
  "92wk_F5" = "slide118_"
) 

# changing names of each element in annotation_list to its corresponding sample id 

names(annotation_list) <- names(prefix_map)


# iterating through annotation_list and adding the proper suffix to the cell ids

annotation_list <- imap(annotation_list, function(df, name) {
  
  df <- annotation_list[[name]] # extracting each dataframe from annotation_list
  
  # appending suffix to cell ids
  
  df <- df %>% 
    mutate(Cell.ID = paste0(prefix_map[[name]], Cell.ID), 
           
           # adding sample id column
           Sample.ID = name)  
  
  return(df)
})


# merging together XE annotations that are on the same slide 
# run 1 slide 0021999
ann_0021999 <- annotation_list[["92wk_M1"]] 

# SANITY CHECK
unique(ann_0021999$Sample.ID)


# run 1 slide 0022573
ann_0022573 <- bind_rows(
  annotation_list[["59wk_F2"]],
  annotation_list[["59wk_M1"]],
  annotation_list[["92wk_F2"]],
  annotation_list[["92wk_M2"]]
) 

# SANITY CHECK
unique(ann_0022573$Sample.ID)



# run 2 slide 0021991
ann_0021991 <- bind_rows(
  annotation_list[["16wk_M1"]], 
  annotation_list[["36wk_F1"]], 
  annotation_list[["36wk_M2"]]  
)

# SANITY CHECK
unique(ann_0021991$Sample.ID)


# run 2 slide 0021998
ann_0021998 <- bind_rows(
  annotation_list[["16wk_F1"]], 
  annotation_list[["16wk_M2"]], 
  annotation_list[["36wk_M1"]]  
)

# SANITY CHECK
unique(ann_0021998$Sample.ID)

# run 4 slide 0069039
ann_0069039 <- bind_rows(
  annotation_list[["16wk_M4"]], 
  annotation_list[["59wk_F5"]], 
  annotation_list[["92wk_F3"]], 
  annotation_list[["92wk_M3"]]  
)

# SANITY CHECK
unique(ann_0069039$Sample.ID)


# run 4 slide 0069045
ann_0069045 <- bind_rows(
  annotation_list[["16wk_F4"]], 
  annotation_list[["36wk_F3"]], 
  annotation_list[["59wk_M3"]]  
)

# SANITY CHECK
unique(ann_0069045$Sample.ID)


# run 5 slide 0069118
ann_0069118 <- bind_rows(
  annotation_list[["16wk_F2"]], 
  annotation_list[["36wk_F4"]], 
  annotation_list[["36wk_M4"]], 
  annotation_list[["59wk_F3"]], 
  annotation_list[["59wk_M5"]], 
  annotation_list[["92wk_F5"]]  
)

# SANITY CHECK 
unique(ann_0069118$Sample.ID)

# set rownames to cell ID so that they match seurat obj cell names 

rownames(ann_0021999) <- ann_0021999$Cell.ID
rownames(ann_0022573) <- ann_0022573$Cell.ID
rownames(ann_0021991) <- ann_0021991$Cell.ID
rownames(ann_0021998) <- ann_0021998$Cell.ID
rownames(ann_0069039) <- ann_0069039$Cell.ID
rownames(ann_0069045) <- ann_0069045$Cell.ID
rownames(ann_0069118) <- ann_0069118$Cell.ID


# reordering rows to match seurat object

ann_0021999 <- ann_0021999[Cells(xenium.obj1), , drop = T]
ann_0022573 <- ann_0022573[Cells(xenium.obj2), , drop = T]
ann_0021991 <- ann_0021991[Cells(xenium.obj3), , drop = T]
ann_0021998 <- ann_0021998[Cells(xenium.obj4), , drop = T]
ann_0069039 <- ann_0069039[Cells(xenium.obj5), , drop = T]
ann_0069045 <- ann_0069045[Cells(xenium.obj6), , drop = T]
ann_0069118 <- ann_0069118[Cells(xenium.obj7), , drop = T]





## 03c. Adding cell ID and sample ID to seurat metadata ----

xenium.obj1 <- AddMetaData(xenium.obj1, ann_0021999$Sample.ID, col.name = "sample_ID")

dim(xenium.obj1)
# 480 374883

xenium.obj2 <- AddMetaData(xenium.obj2, ann_0022573$Sample.ID, col.name = "sample_ID")

dim(xenium.obj2)
# 480 418997

xenium.obj3 <- AddMetaData(xenium.obj3, ann_0021991$Sample.ID, col.name = "sample_ID")

dim(xenium.obj3)
# 480 413608

xenium.obj4 <- AddMetaData(xenium.obj4, ann_0021998$Sample.ID, col.name = "sample_ID")

dim(xenium.obj4)
# 480 380795

xenium.obj5 <- AddMetaData(xenium.obj5, ann_0069039$Sample.ID, col.name = "sample_ID")

dim(xenium.obj5)
# 480 402723

xenium.obj6 <- AddMetaData(xenium.obj6, ann_0069045$Sample.ID, col.name = "sample_ID")

dim(xenium.obj6)
# 480 409647

xenium.obj7 <- AddMetaData(xenium.obj7, ann_0069118$Sample.ID, col.name = "sample_ID")

dim(xenium.obj7)
# 480 415056


# 04. Subset out WT brains in each object ----

## 04a. pull out cells that do not have an NA in their sample_ID metadata
cellskeep_999 <- colnames(xenium.obj1)[!is.na(xenium.obj1$sample_ID)]
cellskeep_573 <- colnames(xenium.obj2)[!is.na(xenium.obj2$sample_ID)]
cellskeep_991 <- colnames(xenium.obj3)[!is.na(xenium.obj3$sample_ID)]
cellskeep_998 <- colnames(xenium.obj4)[!is.na(xenium.obj4$sample_ID)]
cellskeep_039 <- colnames(xenium.obj5)[!is.na(xenium.obj5$sample_ID)]
cellskeep_045 <- colnames(xenium.obj6)[!is.na(xenium.obj6$sample_ID)]
cellskeep_118 <- colnames(xenium.obj7)[!is.na(xenium.obj7$sample_ID)]




## 04b. subset the objects to only WT brains ----
xenium.objWT_999 <- subset(xenium.obj1, cells = cellskeep_999)

dim(xenium.objWT_999)
# 480 54881

xenium.objWT_573 <- subset(xenium.obj2, cells = cellskeep_573)

dim(xenium.objWT_573)
# 480 293814

xenium.objWT_991 <- subset(xenium.obj3, cells = cellskeep_991)

dim(xenium.objWT_991)
# 480 226647

xenium.objWT_998 <- subset(xenium.obj4, cells = cellskeep_998)

dim(xenium.objWT_998)
# 480 205683

xenium.objWT_039 <- subset(xenium.obj5, cells = cellskeep_039)
# 480 251739

xenium.objWT_045 <- subset(xenium.obj6, cells = cellskeep_045)

dim(xenium.objWT_045)
# 480 196482

xenium.objWT_118 <- subset(xenium.obj7, cells = cellskeep_118)

dim(xenium.objWT_118)
# 480 413476

# 05. Save each file with only the WT brains ----

qs_save(xenium.objWT_999, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0021999.qs2")
qs_save(xenium.objWT_573, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0022573.qs2")
qs_save(xenium.objWT_991, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0021991.qs2")
qs_save(xenium.objWT_998, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0021998.qs2")
qs_save(xenium.objWT_039, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0069039.qs2")
qs_save(xenium.objWT_045, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0069045.qs2")
qs_save(xenium.objWT_118, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0069118.qs2")


# 06. Implementing BPCells-backed matrix to each obj before merging ----

## 06a. read in individ objs ----
xenium.objWT_999 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0021999.qs2")
xenium.objWT_573 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0022573.qs2")
xenium.objWT_991 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0021991.qs2")
xenium.objWT_998 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0021998.qs2")
xenium.objWT_039 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0069039.qs2")
xenium.objWT_045 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0069045.qs2")
xenium.objWT_118 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0069118.qs2")

## 06b. write matrix directory for each obj to store the counts ----
write_matrix_dir(
  mat = convert_matrix_type( # ensure the counts are stored as integers and not doubles (aka decimals) to compress efficiently
    LayerData(xenium.objWT_999, layer = "counts"), # RAW COUNTS (not normalized or scaled)
    type = "uint32_t"), # unsigned 32-bit integer
  dir = "./bpcells/slide0021999_counts"
)

write_matrix_dir(
  mat = convert_matrix_type(
    LayerData(xenium.objWT_573, layer = "counts"),
    type = "uint32_t"),
  dir = "./bpcells/slide0022573_counts"
)

write_matrix_dir(
  mat = convert_matrix_type(
    LayerData(xenium.objWT_991, layer = "counts"),
    type = "uint32_t"),
  dir = "./bpcells/slide0021991_counts"
)

write_matrix_dir(
  mat = convert_matrix_type(
    LayerData(xenium.objWT_998, layer = "counts"),
    type = "uint32_t"),
  dir = "./bpcells/slide0021998_counts"
)

write_matrix_dir(
  mat = convert_matrix_type(
    LayerData(xenium.objWT_039, layer = "counts"),
    type = "uint32_t"),
  dir = "./bpcells/slide0069039_counts"
)

write_matrix_dir(
  mat = convert_matrix_type(
    LayerData(xenium.objWT_045, layer = "counts"),
    type = "uint32_t"),
  dir = "./bpcells/slide0069045_counts"
)

write_matrix_dir(
  mat = convert_matrix_type(
    LayerData(xenium.objWT_118, layer = "counts"),
    type = "uint32_t"),
  dir = "./bpcells/slide0069118_counts"
)

## 06c. reload as BPCells-backed matrix and reassign counts assay ----

# open BPCells directory 
bp_counts999 <- open_matrix_dir("./bpcells/slide0021999_counts")

# assign BPCells-backed matrix to obj
xenium.objWT_999[["Xenium"]] <- CreateAssay5Object(counts = bp_counts999)

# repeat with other objs
bp_counts573 <- open_matrix_dir("./bpcells/slide0022573_counts")
xenium.objWT_573[["Xenium"]] <- CreateAssay5Object(counts = bp_counts573)

bp_counts991 <- open_matrix_dir("./bpcells/slide0021991_counts")
xenium.objWT_991[["Xenium"]] <- CreateAssay5Object(counts = bp_counts991)

bp_counts998 <- open_matrix_dir("./bpcells/slide0021998_counts")
xenium.objWT_998[["Xenium"]] <- CreateAssay5Object(counts = bp_counts998)

bp_counts039 <- open_matrix_dir("./bpcells/slide0069039_counts")
xenium.objWT_039[["Xenium"]] <- CreateAssay5Object(counts = bp_counts039)

bp_counts045 <- open_matrix_dir("./bpcells/slide0069045_counts")
xenium.objWT_045[["Xenium"]] <- CreateAssay5Object(counts = bp_counts045)

bp_counts118 <- open_matrix_dir("./bpcells/slide0069118_counts")
xenium.objWT_118[["Xenium"]] <- CreateAssay5Object(counts = bp_counts118)


# ATTACHING UPDATED XE ANNOTATIONS ----
# 4/20/26: some XE full brain annotations contained cells from contralateral hemi (samples 16wk_M1 on slide 0021991, 36wk_F4, 59wk_F3 and 59wk_M5 on slide 0069118)
# replacing sample ID metadata with new XE annotations

xenium.objWT_991 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/Old_objs/20260416_WT_slide0021991.qs2")
xenium.objWT_118 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/Old_objs/20260416_WT_slide0069118.qs2")

# remove old sample_ID metadata col
xenium.objWT_991$sample_ID <- NULL
xenium.objWT_118$sample_ID <- NULL


# paths to the XE full brain annotations (all annotations belonging to this slide including the updated annotations)
files <- list(
  # run 2
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/updated_run2_slide0021991_16wk_M1_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run2_slide0021991_36wk_F1_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run2_slide0021991_36wk_M2_cells_stats.csv",

  # run 5
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run5_slide0069118_16wk_F2_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/updated_run5_slide0069118_36wk_F4_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run5_slide0069118_36wk_M4_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/updated_run5_slide0069118_59wk_F3_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/updated_run5_slide0069118_59wk_M5_cells_stats.csv",
  "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/XE_annotations/fullbrain/run5_slide0069118_92wk_F5_cells_stats.csv"
  
)


# formatting XE full brain annotations to add to seurat metadata 

# reading and tagging each csv cell stats file with corresponding sample name 

annotation_list <- map(files, function(x) { 
  
  # reading each csv file, skipping first 2 metadata lines
  
  read.csv(x, skip = 2, blank.lines.skip = TRUE)
  
})


# defining which samples are in each xenium object for proper prefix placement

prefix_map <- c(
  "16wk_M1" = "slide991_", 
  "36wk_F1" = "slide991_", 
  "36wk_M2" = "slide991_", 
 
  "16wk_F2" = "slide118_", 
  "36wk_F4" = "slide118_",
  "36wk_M4" = "slide118_",
  "59wk_F3" = "slide118_",
  "59wk_M5" = "slide118_",
  "92wk_F5" = "slide118_"
) 

# changing names of each element in annotation_list to its corresponding sample id 

names(annotation_list) <- names(prefix_map)


# iterating through annotation_list and adding the proper suffix to the cell ids

annotation_list <- imap(annotation_list, function(df, name) {
  
  df <- annotation_list[[name]] # extracting each dataframe from annotation_list
  
  # appending suffix to cell ids
  
  df <- df %>% 
    mutate(Cell.ID = paste0(prefix_map[[name]], Cell.ID), 
           
           # adding sample id column
           Sample.ID = name)  
  
  return(df)
})


# merging together XE annotations that are on the same slide 
# run 2 slide 0021991
ann_0021991 <- bind_rows(
  annotation_list[["16wk_M1"]], 
  annotation_list[["36wk_F1"]], 
  annotation_list[["36wk_M2"]]  
)

# SANITY CHECK
unique(ann_0021991$Sample.ID)


# run 5 slide 0069118
ann_0069118 <- bind_rows(
  annotation_list[["16wk_F2"]], 
  annotation_list[["36wk_F4"]], 
  annotation_list[["36wk_M4"]], 
  annotation_list[["59wk_F3"]], 
  annotation_list[["59wk_M5"]], 
  annotation_list[["92wk_F5"]]  
)

# SANITY CHECK 
unique(ann_0069118$Sample.ID)

# set rownames to cell ID so that they match seurat obj cell names 
rownames(ann_0021991) <- ann_0021991$Cell.ID
rownames(ann_0069118) <- ann_0069118$Cell.ID


# reordering rows to match seurat object
ann_0021991 <- ann_0021991[Cells(xenium.objWT_991), , drop = T]
ann_0069118 <- ann_0069118[Cells(xenium.objWT_118), , drop = T]

# Adding cell ID and sample ID to seurat metadata 

xenium.objWT_991 <- AddMetaData(xenium.objWT_991, ann_0021991$Sample.ID, col.name = "sample_ID")

dim(xenium.objWT_991)
# 480 226647

xenium.objWT_118 <- AddMetaData(xenium.objWT_118, ann_0069118$Sample.ID, col.name = "sample_ID")

dim(xenium.objWT_118)
# 480 413476


## 04a. pull out cells that do not have an NA in their sample_ID metadata
cellskeep_991 <- colnames(xenium.objWT_991)[!is.na(xenium.objWT_991$sample_ID)]
cellskeep_118 <- colnames(xenium.objWT_118)[!is.na(xenium.objWT_118$sample_ID)]


## subset the objects to only WT brains 
xenium.objWT_991 <- subset(xenium.objWT_991, cells = cellskeep_991)

dim(xenium.objWT_991)
# 480 215712

xenium.objWT_118 <- subset(xenium.objWT_118, cells = cellskeep_118)

dim(xenium.objWT_118)

# 480 400924


# delete old directory 
unlink("./bpcells/slide0021991_counts", recursive = TRUE)
unlink("./bpcells/slide0069118_counts", recursive = TRUE)



# Save updated objs 

qs_save(xenium.objWT_991, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260420_WT_slide0021991.qs2")
qs_save(xenium.objWT_118, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260420_WT_slide0069118.qs2")





























# 07. Merge objects and saving merged file ----

## read in objs ----

xenium.objWT_999 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0021999.qs2")
xenium.objWT_573 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0022573.qs2")
xenium.objWT_991 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260420_WT_slide0021991.qs2")
xenium.objWT_998 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0021998.qs2")
xenium.objWT_039 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0069039.qs2")
xenium.objWT_045 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260416_WT_slide0069045.qs2")
xenium.objWT_118 <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/WTsubset_individ_seuratobjs/20260420_WT_slide0069118.qs2")




## 07a. make list of all individual WT subsetted objs ----
xenium_list <- c(xenium.objWT_999, 
                 xenium.objWT_573, 
                 xenium.objWT_991, 
                 xenium.objWT_998, 
                 xenium.objWT_039, 
                 xenium.objWT_045, 
                 xenium.objWT_118)

## 07b. merge objets using the list created above ----
WTmerged.obj <- merge(xenium_list[[1]], y = xenium_list[-1])

# SANITY CHECK
dim(WTmerged.obj)
# 480 1619235

## 07c. Join layers ----
WTmerged.obj <- JoinLayers(WTmerged.obj, add.prefix = FALSE)

## 07d. subset nCount > 5 ----
WTmerged.obj <- subset(WTmerged.obj, subset = nCount_Xenium > 5)

# SANITY CHECK
dim(WTmerged.obj)
# 480 1609408

# 08. save merged object ----

qs_save(WTmerged.obj, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260420_mergedobj_no_metadata.qs2")

# 09. Attaching metadata from xlsx file ----

## 09a. load merged object ----

WTmerged.obj <- qs_read("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260420_mergedobj_no_metadata.qs2")

## 09b. read in metadata excel sheet ----
WT_metadata <- data.frame(
  read_excel("/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/metadata/WT_AgingXMetabolism_Xenium_metadata.xlsx"))

## 09c. extract metadata from merged seurat obj ----
meta <- WTmerged.obj@meta.data

## 09d. convert rownames to cell id ----
meta <- meta %>%
  rownames_to_column(var = "cell_id")

## 09e. join metadata excel sheet to merged obj metadata 
meta <- meta %>%
  left_join(WT_metadata, by = "sample_ID"
            #, keep = T
  )

## 09f. revert cell id back to columns ----
meta <- column_to_rownames(meta, var = "cell_id")

## 09g. attach metadata back to merged obj ----
WTmerged.obj@meta.data <- meta


# SANITY CHECK
dim(WTmerged.obj)

table(WTmerged.obj$sample_ID, WTmerged.obj$slide)

table(WTmerged.obj$Age, WTmerged.obj$Sex)

table(WTmerged.obj$sample_ID, WTmerged.obj$Age)

# 10. save final part 1 merged obj ----
qs_save(WTmerged.obj, "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/20260417_mergedobjv1_01.qs2")
