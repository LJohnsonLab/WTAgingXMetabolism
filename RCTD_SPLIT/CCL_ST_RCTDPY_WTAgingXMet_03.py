'''
Name: CCL_ST_RCTDPY_WTAgingXMet_03.py
Project: Wildtype Aging Cerebral Metabolism
Purpose: Analysis of 4 different age groups (16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data  
        -  03.  running RCTD-py  
Input Files:   
Final Output Files: 
Kernel: scverse_env (Python 3.11.14)
Date created: 8/5/26
Last updated: 8/5/26
Author: Chloe Lucido
'''

from rctd import Reference, run_rctd
from rctd._types import RCTDConfig
import anndata as ad
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import spatialdata as sd
from spatialdata_io import xenium
from pathlib import Path
import h5py
from scipy.sparse import issparse

# settings for saving as h5ad 
pd.set_option("future.infer_string", False) # prevents pandas from making pyarrow-backed string columns in first place 
pd.options.mode.string_storage = "python"
ad.settings.allow_write_nullable_strings = True

###### PATHS ######
spatial_obj = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/RCTD_SPLIT/20260722_merged.h5ad"
ref_obj = "/Users/cclu223/Desktop/ABC_reference/WMB_10xv3_data/downsampled_objs/20260723_raw_WMB_10xv3_FINAL.h5ad"
output_path = "/Users/cclu223/Desktop/WT_AgingXMet/Xenium_Analysis/obj_files/RCTD_SPLIT/"
DATE = "20260805_" 

export_save_dir = Path(output_path) / f"{DATE}rctd_results_xe_Subcluster"
###################

# 01. read objs
spatial = ad.read_h5ad(spatial_obj)
reference = Reference(ad.read_h5ad(ref_obj), cell_type_col = "Subcluster")

# 02. set class_df configuration

## 02a. create class_df
class_df = {
    "NK Cells": "NK Cells",
    "T Cells": "T Cells",
    "Microglia": "Microglia",
    "BAM": "BAM",
    "Monocytes": "Monocytes",
    "DC": "DC",
    "ABC": "ABC",
    "VLMC": "VLMC",
    "Pericytes": "Pericytes",
    "SMC": "SMC",
    "Endothelial Cells": "Endothelial Cells",
    "OPC": "OPC",
    "Oligodendrocytes": "Oligodendrocytes",
    "Astrocytes": "Astrocytes",
    "Astroependymal": "Astroependymal",
    "Tanycytes": "Tanycytes",
    "Ependymal": "Ependymal",
    "Hypendymal": "Hypendymal",
    "CP": "CP",
    "Dopaminergic_Neurons": "Dopaminergic_Neurons",
    "Serotonergic_Neurons": "Serotonergic_Neurons",
    "LAMP5": "GABAergic Neuron",
    "PVALB": "GABAergic Neuron",
    "SNCG": "GABAergic Neuron",
    "SST": "GABAergic Neuron",
    "VIP": "GABAergic Neuron",
    "L2_3_IT": "Glutamatergic Neuron",
    "L5_IT": "Glutamatergic Neuron",
    "L6_IT": "Glutamatergic Neuron",
    "L5_ET": "Glutamatergic Neuron",
    "L6_CT": "Glutamatergic Neuron"
}

config = RCTDConfig(class_df = class_df)
result = run_rctd(spatial, reference, mode = "doublet", config = config)



###### 03. Export result for SPLIT's reconstruct_rctd_from_rctdpy() (R) ######
#
# from Claude (8/5/26)
#
# NOTE: this assumes `result` exposes the attribute names documented in
# reconstruct_rctd_from_rctdpy()'s own docstring (weights, weights_doublet,
# spot_class, first_type, second_type, first_class, second_class, min_score,
# singlet_score, pixel_mask, cell_type_names). I haven't been able to verify
# these against your actual installed `rctd` package version -- if this
# throws an AttributeError, run `vars(result)` or `dir(result)` and tell me
# what's actually there so I can fix the mapping.
#
# cell_type_col="Subcluster" matches what you confirmed earlier as the
# column in your reference h5ad holding the 31 fine-grained cell type labels
# (the same ones used as keys in class_df above).
 
def export_rctdpy_for_reconstruction(xe, result, reference_h5ad, cell_type_col, save_dir):
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
 
    n_weights = result.weights.shape[0]
    if hasattr(result, "pixel_mask"):
        n_mask_true = int(np.sum(result.pixel_mask))
        if n_mask_true != n_weights:
            raise ValueError(
                f"result.pixel_mask has {n_mask_true} True entries but "
                f"result.weights has {n_weights} rows -- these must match."
            )
        cell_ids = xe.obs_names[result.pixel_mask]
    else:
        print("WARNING: result has no `pixel_mask` attribute. Assuming "
              "weights cover all of xe.obs_names in order -- verify this "
              "is actually true before trusting the output.")
        cell_ids = xe.obs_names
        if len(cell_ids) != n_weights:
            raise ValueError(
                f"len(xe.obs_names)={len(cell_ids)} != "
                f"result.weights.shape[0]={n_weights}; cannot safely assume "
                f"no filtering was applied."
            )
 
    assert len(cell_ids) == n_weights
    assert len(cell_ids) == len(result.spot_class)
 
    cell_type_names = list(result.cell_type_names)
 
    pd.DataFrame({"cell_id": cell_ids}).to_parquet(
        save_dir / "cell_ids.parquet", index=False
    )
 
    pd.DataFrame(
        np.asarray(result.weights), columns=cell_type_names
    ).to_parquet(save_dir / "weights.parquet")
 
    pd.DataFrame(
        np.asarray(result.weights_doublet), columns=["w_1", "w_2"]
    ).to_parquet(save_dir / "weights_doublet.parquet")
 
    def resolve_type_names(arr):
        arr = np.asarray(arr)
        if np.issubdtype(arr.dtype, np.integer):
            return [cell_type_names[i] for i in arr]
        return list(arr)
 
    spot_df = pd.DataFrame({
        "cell_id":          cell_ids,
        "spot_class":       result.spot_class,
        "first_type":       result.first_type,
        "second_type":      result.second_type,
        "first_class":      result.first_class,
        "second_class":     result.second_class,
        "min_score":        result.min_score,
        "singlet_score":    result.singlet_score,
        "first_type_name":  resolve_type_names(result.first_type),
        "second_type_name": resolve_type_names(result.second_type),
    })
    spot_df.to_parquet(save_dir / "spot_results.parquet")
 
    pixel_mask_full = (
        np.asarray(result.pixel_mask) if hasattr(result, "pixel_mask")
        else np.ones(len(xe.obs_names), dtype=bool)
    )
    pd.DataFrame({
        "cell_id": xe.obs_names,
        "pixel_mask": pixel_mask_full,
    }).to_parquet(save_dir / "pixel_mask.parquet")
 
    pd.DataFrame({"cell_type_names": cell_type_names}).to_parquet(
        save_dir / "metadata.parquet"
    )
 
    # reference_profiles.h5 -- built directly from your existing reference
    # h5ad rather than a live rctd-py Reference object. This is a plain
    # per-cell-type mean over all reference cells of that type; if rctd-py's
    # internal fitting applied cell_min / n_max_cells filtering on the
    # reference before building its own profiles, this won't be
    # bit-identical -- shouldn't meaningfully change SPLIT's output, but
    # worth knowing.
    ref = ad.read_h5ad(reference_h5ad)
    if cell_type_col not in ref.obs.columns:
        raise ValueError(
            f"'{cell_type_col}' not in reference obs columns: "
            f"{list(ref.obs.columns)}"
        )
 
    X = ref.X
    if issparse(X):
        X = X.toarray()
    X = pd.DataFrame(X, index=ref.obs_names, columns=ref.var_names)
    profiles = X.groupby(ref.obs[cell_type_col].values).mean()  # cell_types x genes
 
    missing = set(cell_type_names) - set(profiles.index)
    if missing:
        raise ValueError(
            f"Cell types in rctd-py results but not found in reference's "
            f"'{cell_type_col}' column: {missing}"
        )
    profiles = profiles.loc[cell_type_names]
 
    with h5py.File(save_dir / "reference_profiles.h5", "w") as f:
        f.create_dataset("profiles", data=profiles.to_numpy())
        f.create_dataset(
            "cell_type_names",
            data=np.asarray(profiles.index.tolist(), dtype="S"),
        )
        f.create_dataset(
            "gene_names", data=np.asarray(profiles.columns.tolist(), dtype="S")
        )
 
    print(f"Wrote reconstruct_rctd_from_rctdpy()-compatible files to {save_dir} "
          f"({n_weights} cells, {len(cell_type_names)} cell types)")
 
 
export_rctdpy_for_reconstruction(
    xe=spatial,
    result=result,
    reference_h5ad=ref_obj,
    cell_type_col="Subcluster",
    save_dir=str(export_save_dir),
)