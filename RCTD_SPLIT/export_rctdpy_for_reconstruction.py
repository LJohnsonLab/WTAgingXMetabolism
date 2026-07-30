"""
Export rctd-py script/library-API results into the file layout expected by
SPLIT's reconstruct_rctd_from_rctdpy() R function.

This is a direct implementation of the Python snippet documented in that
function's own docstring -- use this instead of the CLI's h5ad output when
you ran rctd-py via its library API directly (which doesn't autosave).

IMPORTANT -- verify before trusting this blindly:
  This assumes your `result` object exposes the same attribute names as the
  docstring's example (`weights`, `weights_doublet`, `spot_class`,
  `first_type`, `second_type`, `first_class`, `second_class`, `min_score`,
  `singlet_score`, `pixel_mask`, `cell_type_names`), and that your
  `reference` object exposes `profiles`, `cell_type_names`, `gene_names`.
  Run `vars(result)` / `dir(result)` and `vars(reference)` / `dir(reference)`
  first to confirm these exist under these exact names in your installed
  rctd-py version -- attribute names can drift between versions and I can't
  verify this against your actual objects.

  This also assumes doublet mode specifically (reconstruct_rctd_from_rctdpy
  is doublet-mode only, matching SPLIT's own doublet-mode requirement).
"""

import h5py
import numpy as np
import pandas as pd
from pathlib import Path


def export_rctdpy_for_reconstruction(
    xe,                      # AnnData: your spatial object (all original cells)
    result,                  # rctd-py DoubletResult-like object from run_rctd()
    reference_h5ad: str,     # path to your existing reference scRNA-seq h5ad
    cell_type_col: str,      # obs column in reference_h5ad with cell type labels
    save_dir: str,
):
    """
    Writes the 6 files reconstruct_rctd_from_rctdpy() expects into save_dir:
    cell_ids.parquet, weights.parquet, weights_doublet.parquet,
    spot_results.parquet, pixel_mask.parquet, metadata.parquet,
    reference_profiles.h5.

    Reference profiles are built directly from your existing reference.h5ad
    (per-cell-type mean expression), rather than pulled from a live rctd-py
    Reference object. NOTE: this is a plain per-type mean over all reference
    cells of that type -- if rctd-py's internal fitting applied cell_min /
    n_max_cells filtering on the reference before building its own profiles,
    this won't be bit-identical to what RCTD used internally. Shouldn't
    meaningfully change SPLIT's output, but worth knowing.
    """
    import anndata as ad
    from scipy.sparse import issparse

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    # --- sanity checks up front, so failures are legible rather than cryptic ---
    n_weights = result.weights.shape[0]
    if hasattr(result, "pixel_mask"):
        n_mask_true = int(np.sum(result.pixel_mask))
        if n_mask_true != n_weights:
            raise ValueError(
                f"result.pixel_mask has {n_mask_true} True entries but "
                f"result.weights has {n_weights} rows -- these must match. "
                f"Check whether pixel_mask is boolean over ALL xe.obs_names "
                f"(expected) vs already subset."
            )
        cell_ids = xe.obs_names[result.pixel_mask]
    else:
        # Fall back: assume no filtering was applied, i.e. every cell in xe
        # is present in result.weights, in the same order.
        print("WARNING: result has no `pixel_mask` attribute. Assuming "
              "weights cover all of xe.obs_names in order -- verify this "
              "is actually true for your rctd-py version before proceeding.")
        cell_ids = xe.obs_names
        if len(cell_ids) != n_weights:
            raise ValueError(
                f"len(xe.obs_names)={len(cell_ids)} != "
                f"result.weights.shape[0]={n_weights}; cannot safely assume "
                f"no filtering was applied. You need the actual pixel mask."
            )

    assert len(cell_ids) == n_weights
    assert len(cell_ids) == len(result.spot_class)

    cell_type_names = list(result.cell_type_names)

    # --- cell_ids.parquet ---
    pd.DataFrame({"cell_id": cell_ids}).to_parquet(
        save_dir / "cell_ids.parquet", index=False
    )

    # --- weights.parquet (cells x cell_types) ---
    pd.DataFrame(
        np.asarray(result.weights), columns=cell_type_names
    ).to_parquet(save_dir / "weights.parquet")

    # --- weights_doublet.parquet (cells x 2, columns MUST be w_1/w_2) ---
    pd.DataFrame(
        np.asarray(result.weights_doublet), columns=["w_1", "w_2"]
    ).to_parquet(save_dir / "weights_doublet.parquet")

    # --- spot_results.parquet ---
    # first_type/second_type may be integer indices (need lookup into
    # cell_type_names) or already resolved strings, depending on your
    # rctd-py version. Handle both.
    def resolve_type_names(arr):
        arr = np.asarray(arr)
        if np.issubdtype(arr.dtype, np.integer):
            return [cell_type_names[i] for i in arr]
        return list(arr)  # already strings

    df = pd.DataFrame({
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
    df.to_parquet(save_dir / "spot_results.parquet")

    # --- pixel_mask.parquet (over ALL original obs_names) ---
    pixel_mask_full = (
        np.asarray(result.pixel_mask) if hasattr(result, "pixel_mask")
        else np.ones(len(xe.obs_names), dtype=bool)
    )
    pd.DataFrame({
        "cell_id": xe.obs_names,
        "pixel_mask": pixel_mask_full,
    }).to_parquet(save_dir / "pixel_mask.parquet")

    # --- metadata.parquet ---
    pd.DataFrame({"cell_type_names": cell_type_names}).to_parquet(
        save_dir / "metadata.parquet"
    )

    # --- reference_profiles.h5 (built from your existing reference.h5ad) ---
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
    profiles = profiles.loc[cell_type_names]  # reorder to match weights columns

    with h5py.File(save_dir / "reference_profiles.h5", "w") as f:
        f.create_dataset("profiles", data=profiles.to_numpy())  # cell_types x genes
        f.create_dataset(
            "cell_type_names",
            data=np.asarray(profiles.index.tolist(), dtype="S"),
        )
        f.create_dataset(
            "gene_names", data=np.asarray(profiles.columns.tolist(), dtype="S")
        )

    print(f"Wrote reconstruct_rctd_from_rctdpy()-compatible files to {save_dir} "
          f"({n_weights} cells, {len(cell_type_names)} cell types)")


if __name__ == "__main__":
    # Example usage -- replace with your actual objects
    # export_rctdpy_for_reconstruction(
    #     xe=xe,
    #     result=xe_rctd_results,
    #     reference_h5ad="reference.h5ad",
    #     cell_type_col="cell_type",
    #     save_dir="rctd_results_xe_level2_chromium",
    # )
    pass
