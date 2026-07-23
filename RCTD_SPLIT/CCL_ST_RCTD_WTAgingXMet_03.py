########################################################################
# Name: CCL_ST_RCTD_WTAgingXMet_03.py
# Project: Wildtype Aging Cerebral Metabolism
# Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
#         -  03.  running RCTD 
# Input Files:   
# Final Output Files: 
# Kernel: scverse_env (Python 3.11.14)
# Date created: 7/22/26
# Last updated: 7/23/26
# Author: Chloe Lucido
########################################################################

from rctd import Reference, run_rctd
import anndata as ad
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import spatialdata as sd
from spatialdata_io import xenium

# settings for saving as h5ad 
pd.set_option("future.infer_string", False)
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True


###### PATHS ######
ref_obj_path = "/Users/cclu223/Desktop/ABC_reference/WMB_10xv3_data/downsampled_objs/20260714_log2_WMB_10xv3_FINAL.h5ad"
spatial_obj_path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/RCTD_SPLIT/20260722_merged.h5ad"
# = ""
###################

reference = ad.read_h5ad(ref_obj_path)

print(reference.obs.columns)
print(reference.obs["Subcluster"].isna().sum())
print(reference.obs["Broad_Cluster"].isna().sum())


# RCTD

reference = Reference(reference, cell_type_col="Subcluster") # CHECK CELL TYPE COLUMN !!!

spatial = ad.read_h5ad(spatial_obj_path)

# Run RCTD — handles normalization, sigma estimation, and deconvolution
result = run_rctd(spatial, reference, mode="doublet")