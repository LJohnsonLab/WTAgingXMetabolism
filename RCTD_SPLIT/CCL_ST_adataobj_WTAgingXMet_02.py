'''
Name: CCL_ST_adataobj_WTAgingXMet_02.py
Project: Wildtype Aging Cerebral Metabolism
Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
        -  02.  building anndata obj for RCTD analysis   
Input Files:   
Final Output Files: 
Kernel: scverse_env (Python 3.11.14)
Date created: 7/16/26
Last updated: 7/28/26
Author: Chloe Lucido
'''

from rctd import Reference, run_rctd
import anndata as ad
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import spatialdata as sd
from spatialdata_io import xenium
from pathlib import Path

# settings for saving as h5ad 
pd.set_option("future.infer_string", False) # prevents pandas from making pyarrow-backed string columns in first place 
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

###### PATHS ######
slide999_path = '/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run1_20250516__201914__20250516_AgingXMetabolism_1/Slide1_output-XETG00118__0021999__Region_1__20250516__201923'
slide573_path = '/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run1_20250516__201914__20250516_AgingXMetabolism_1/Slide2_output-XETG00118__0022573__Region_1__20250516__201924'
slide991_path = '/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run2_20250606__202944__20250606_AgingXMetabolism_2/Slide1_output-XETG00118__0021991__Region_1__20250606__202953'
slide998_path = '/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run2_20250606__202944__20250606_AgingXMetabolism_2/Slide2_output-XETG00118__0021998__Region_1__20250606__202953'
slide039_path = '/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run4_20250919__205218__20250919_AgingXMetabolism_4/Slide1_output-XETG00118__0069039__Region_1__20250919__205228'
slide045_path = '/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run4_20250919__205218__20250919_AgingXMetabolism_4/Slide2_output-XETG00118__0069045__Region_1__20250919__205229'
slide118_path = '/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run5_20260306_210910_20260306_AgingXMetabolism_5/slide1_output-XETG00118__0069118__Region_1__20260306__210910'
seurat_metadata = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/csv/metadata/20260722_mergedobj_metadata_for_pyRCTD.csv"
output_path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/RCTD_SPLIT/"
DATE = "20260722_" 

merged_adata_path = Path(OUTPUT_DIR) / f"{DATE}merged.h5ad"
###################

# 01. build WTAgingXMet anndata obj ----
## 1a. create sdata objs for each slide ----
sdata999 = xenium(slide999_path, cells_as_circles = True)
sdata573 = xenium(slide573_path, cells_as_circles = True)
sdata991 = xenium(slide991_path, cells_as_circles = True)
sdata998 = xenium(slide998_path, cells_as_circles = True)
sdata039 = xenium(slide039_path, cells_as_circles = True)
sdata045 = xenium(slide045_path, cells_as_circles = True)
sdata118 = xenium(slide118_path, cells_as_circles = True)

## 1b.  set index to cell ids and add slide id prefix ################################################
### SLIDE 999--------------

adata999 = sdata999["table"]

#### set index (align)  
adata999.obs_names = [f"slide999_{bc}" for bc in adata999.obs["cell_id"]]

#### contatenate index to cell_id so they are aligned 
adata999.obs["cell_id"] = adata999.obs_names

#### make sure datatype is compatible with seurat metadata
adata999.obs.index = adata999.obs.index.astype(str)

#### SANITY CHECK 
print("Slide prefix added: \n", adata999.obs["cell_id"].head(), "\n")

### SLIDE 573 --------------

adata573 = sdata573["table"]

#### set index (align)  
adata573.obs_names = [f"slide573_{bc}" for bc in adata573.obs["cell_id"]]

#### contatenate index to cell_id so they are aligned 
adata573.obs["cell_id"] = adata573.obs_names

#### make sure datatype is compatible with seurat metadata
adata573.obs.index = adata573.obs.index.astype(str)

#### SANITY CHECK 
print("Slide prefix added: \n", adata573.obs["cell_id"].head(), "\n")


### SLIDE 991 --------------

adata991 = sdata991["table"]

#### set index (align)  
adata991.obs_names = [f"slide991_{bc}" for bc in adata991.obs["cell_id"]]

#### contatenate index to cell_id so they are aligned 
adata991.obs["cell_id"] = adata991.obs_names

#### make sure datatype is compatible with seurat metadata
adata991.obs.index = adata991.obs.index.astype(str)

#### SANITY CHECK 
print("Slide prefix added: \n", adata991.obs["cell_id"].head(), "\n")

### SLIDE 998 ------------

adata998 = sdata998["table"]

#### set index (align)  
adata998.obs_names = [f"slide998_{bc}" for bc in adata998.obs["cell_id"]]

#### contatenate index to cell_id so they are aligned 
adata998.obs["cell_id"] = adata998.obs_names

#### make sure datatype is compatible with seurat metadata
adata998.obs.index = adata998.obs.index.astype(str)

#### SANITY CHECK 
print("Slide prefix added: \n", adata998.obs["cell_id"].head(), "\n")

### SLIDE 039 ------------

adata039 = sdata039["table"]

#### set index (align)  
adata039.obs_names = [f"slide039_{bc}" for bc in adata039.obs["cell_id"]]

#### contatenate index to cell_id so they are aligned 
adata039.obs["cell_id"] = adata039.obs_names

#### make sure datatype is compatible with seurat metadata
adata039.obs.index = adata039.obs.index.astype(str)

#### SANITY CHECK 
print("Slide prefix added: \n", adata039.obs["cell_id"].head(), "\n")

### SLIDE 045 -----------

adata045 = sdata045["table"]

#### set index (align)  
adata045.obs_names = [f"slide045_{bc}" for bc in adata045.obs["cell_id"]]

#### contatenate index to cell_id so they are aligned 
adata045.obs["cell_id"] = adata045.obs_names

#### make sure datatype is compatible with seurat metadata
adata045.obs.index = adata045.obs.index.astype(str)

#### SANITY CHECK 
print("Slide prefix added: \n", adata045.obs["cell_id"].head(), "\n")


### SLIDE 118 -----------

adata118 = sdata118["table"]

#### set index (align)  
adata118.obs_names = [f"slide118_{bc}" for bc in adata118.obs["cell_id"]]

#### contatenate index to cell_id so they are aligned 
adata118.obs["cell_id"] = adata118.obs_names

#### make sure datatype is compatible with seurat metadata
adata118.obs.index = adata118.obs.index.astype(str)

#### SANITY CHECK 
print("Slide prefix added: \n", adata118.obs["cell_id"].head(), "\n")

## 1c.  load and organize seurat metadata ################################################
### read metadata csv file 
meta_py = pd.read_csv(seurat_metadata)

### SANITY CHECK: checking to see if indices match 
print("From Seurat: \n", meta_py.columns, "\n")


### set index of seurat metadata to cell id
meta_py = meta_py.set_index("cell_id")

### make sure datatype of index is compatible
meta_py.index = meta_py.index.astype(str)


### subset seurat metadata by slide 
meta_slide999 = meta_py[meta_py.index.str.startswith("slide999_")]
meta_slide573 = meta_py[meta_py.index.str.startswith("slide573_")]
meta_slide991 = meta_py[meta_py.index.str.startswith("slide991_")]
meta_slide998 = meta_py[meta_py.index.str.startswith("slide998_")]
meta_slide039 = meta_py[meta_py.index.str.startswith("slide039_")]
meta_slide045 = meta_py[meta_py.index.str.startswith("slide045_")]
meta_slide118 = meta_py[meta_py.index.str.startswith("slide118_")]


## 1d. subset out cells with no seurat metadata ################################################
### subset to WTAgingXMet brains using seurat metadata
adata999_subset = adata999[adata999.obs.index.isin(meta_slide999.index)].copy()
adata573_subset = adata573[adata573.obs.index.isin(meta_slide573.index)].copy()
adata991_subset = adata991[adata991.obs.index.isin(meta_slide991.index)].copy()
adata998_subset = adata998[adata998.obs.index.isin(meta_slide998.index)].copy()
adata039_subset = adata039[adata039.obs.index.isin(meta_slide039.index)].copy()
adata045_subset = adata045[adata045.obs.index.isin(meta_slide045.index)].copy()
adata118_subset = adata118[adata118.obs.index.isin(meta_slide118.index)].copy()

# SANITY CHECK: ensure alignment of cell ids btwn seurat metadata and adata objs
print("Number of common cell_ids btwn seurat metadata and python adata999: ", adata999_subset.obs_names.isin(meta_slide999.index).sum())
print("Number of common cell_ids btwn seurat metadata and python adata573: ", adata573_subset.obs_names.isin(meta_slide573.index).sum())
print("Number of common cell_ids btwn seurat metadata and python adata991: ", adata991_subset.obs_names.isin(meta_slide991.index).sum())
print("Number of common cell_ids btwn seurat metadata and python adata998: ", adata998_subset.obs_names.isin(meta_slide998.index).sum())
print("Number of common cell_ids btwn seurat metadata and python adata039: ", adata039_subset.obs_names.isin(meta_slide039.index).sum())
print("Number of common cell_ids btwn seurat metadata and python adata045: ", adata045_subset.obs_names.isin(meta_slide045.index).sum())
print("Number of common cell_ids btwn seurat metadata and python adata118: ", adata118_subset.obs_names.isin(meta_slide118.index).sum())

print("Seurat metadata number of obs: ", len(meta_py.index)) # all obs from adata objs should sum to this


## 1e. attach seurat metadata to adata objs ################################################
### drop overlapping columns to prevent merging errors 
adata999_subset.obs = adata999_subset.obs.drop(columns = meta_slide999.columns.intersection(adata999_subset.obs.columns))
adata573_subset.obs = adata573_subset.obs.drop(columns = meta_slide573.columns.intersection(adata573_subset.obs.columns))
adata991_subset.obs = adata991_subset.obs.drop(columns = meta_slide991.columns.intersection(adata991_subset.obs.columns))
adata998_subset.obs = adata998_subset.obs.drop(columns = meta_slide998.columns.intersection(adata998_subset.obs.columns))
adata039_subset.obs = adata039_subset.obs.drop(columns = meta_slide039.columns.intersection(adata039_subset.obs.columns))
adata045_subset.obs = adata045_subset.obs.drop(columns = meta_slide045.columns.intersection(adata045_subset.obs.columns))
adata118_subset.obs = adata118_subset.obs.drop(columns = meta_slide118.columns.intersection(adata118_subset.obs.columns))


### join metadata safely
adata999_subset.obs = adata999_subset.obs.join(meta_slide999)
adata573_subset.obs = adata573_subset.obs.join(meta_slide573)
adata991_subset.obs = adata991_subset.obs.join(meta_slide991)
adata998_subset.obs = adata998_subset.obs.join(meta_slide998)
adata039_subset.obs = adata039_subset.obs.join(meta_slide039)
adata045_subset.obs = adata045_subset.obs.join(meta_slide045)
adata118_subset.obs = adata118_subset.obs.join(meta_slide118)

### SANITY CHECK
print(adata999_subset.obs.head())

## 1f. merge adata objs together ################################################ 
### create list of adata objs to merge
adatas = [adata999_subset, adata573_subset, adata991_subset, adata998_subset, adata039_subset, adata045_subset, adata118_subset]

### merge all objs together 
merged_adata = ad.concat(adatas, join="outer", axis=0, label="slide_id", keys=["slide999", "slide573", "slide991", "slide998", "slide039", "slide045", "slide118"], index_unique=None)

'''
NOTES:
!!!next time merge all adata objs together first, then subset and attach seurat metadata 
because I can use the index_unique argmument in ad.concat to assign slide prefixes to cell ids !!!
- join argument can be set to "outer" to keep all genes in the merged obj, or "inner" to keep only genes that are shared across all slides
'''

### SANITY CHECK
print("Merged adata obj: \n", merged_adata.n_obs)
print(merged_adata.obs.columns)


## 1h. save merged adata obj ################################################
merged_adata.obs.index = merged_adata.obs.index.astype(str)

# and re-check/convert any lingering pyarrow-backed columns in both obs and var
merged_adata.obs = merged_adata.obs.apply(
    lambda col: col.astype(object) if str(col.dtype).startswith('string') else col
)
merged_adata.var = merged_adata.var.apply(
    lambda col: col.astype(object) if str(col.dtype).startswith('string') else col
)

merged_adata.write_h5ad(merged_adata_path)
