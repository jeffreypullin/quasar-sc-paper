#!/usr/bin/env python3
"""Harmonise Cytoimmgen h5ad for the QuASAR pipeline.

Writes individual, sex (1=male, 2=female), age, batch, timepoint, and
cell_label = "<Compartment> <timepoint>" so existing cell_label_mask /
COMBINED_MEMBERS selection works. Joins pseudotime and division_state from
the CellDivider predictions CSV by barcode.

Usage:
  harmonise-cytoimmgen-h5ad.py <h5ad> <predictions.csv> <out.h5ad> \\
      [individual=col] [compartment=col] [timepoint=col] [batch=col] \\
      [age=col] [sex=col] [pseudotime=col]
"""

from __future__ import annotations

import sys

import pandas as pd
import scanpy as sc
import scipy.sparse as sp

INDIV_ALIASES = ["individual", "donor_id", "donor", "sample_id", "Donor"]
COMPARTMENT_ALIASES = [
    "compartment",
    "cell_type",
    "cell_label",
    "annotation",
    "CellType",
]
TIMEPOINT_ALIASES = ["timepoint", "time", "Timepoint", "Time", "day", "Day"]
BATCH_ALIASES = ["batch", "culture_batch", "Batch", "exp_batch", "lib_batch"]
AGE_ALIASES = ["age", "Age"]
SEX_ALIASES = ["sex", "Sex"]
PSEUDOTIME_ALIASES = ["pseudotime", "Pseudotime", "PT", "dpt_pseudotime"]
SYMBOL_ALIASES = ["GeneSymbol", "feature_name", "gene_symbols", "symbol", "name"]

SEX_MAP = {
    "male": 1,
    "m": 1,
    "1": 1,
    "1.0": 1,
    "female": 2,
    "f": 2,
    "2": 2,
    "2.0": 2,
    "pato:0000384": 1,
    "pato:0000383": 2,
}

COMPARTMENT_MAP = {
    "naive": "Naive",
    "cd4 naive": "Naive",
    "cd4_naive": "Naive",
    "tn": "Naive",
    "memory": "Memory",
    "cd4 memory": "Memory",
    "cd4_memory": "Memory",
    "tcm": "Memory",
    "tem": "Memory",
}


def first_col(frame: pd.DataFrame, names: list[str]) -> str | None:
    for name in names:
        if name in frame.columns:
            return name
    return None


def parse_overrides(args: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for arg in args:
        if "=" not in arg:
            continue
        key, val = arg.split("=", 1)
        out[key.strip()] = val.strip()
    return out


def resolve_col(
    frame: pd.DataFrame, preferred: str | None, aliases: list[str], label: str
) -> str:
    if preferred and preferred in frame.columns:
        return preferred
    found = first_col(frame, aliases)
    if found is None:
        raise SystemExit(f"Could not find {label} column. Tried {aliases}")
    return found


def normalise_sex(series: pd.Series) -> pd.Series:
    mapped = series.astype(str).str.strip().str.lower().map(SEX_MAP)
    if mapped.isna().any():
        bad = series[mapped.isna()].unique()[:5]
        raise SystemExit(f"Unmapped sex values: {bad}")
    return mapped.astype(int)


def normalise_compartment(series: pd.Series) -> pd.Series:
    raw = series.astype(str).str.strip()
    mapped = raw.str.lower().map(COMPARTMENT_MAP)
    # Keep already-canonical Naive / Memory.
    mapped = mapped.fillna(
        raw.where(raw.isin(["Naive", "Memory"]))
    )
    if mapped.isna().any():
        bad = raw[mapped.isna()].unique()[:10]
        raise SystemExit(f"Unmapped compartment values: {bad}")
    return mapped.astype(str)


def normalise_timepoint(series: pd.Series) -> pd.Series:
    raw = series.astype(str).str.strip()
    # 16 -> 16h, 16h -> 16h, day0 -> day0
    out = []
    for val in raw:
        if val.endswith("h") or val.startswith("day") or val.startswith("Day"):
            out.append(val if not val[0].isdigit() or val.endswith("h") else f"{val}h")
        elif val.replace(".", "", 1).isdigit():
            num = str(int(float(val))) if "." not in val else val
            out.append(f"{num}h")
        else:
            out.append(val)
    return pd.Series(out, index=series.index)


def strip_ensembl_version(ids) -> pd.Index:
    return pd.Index([str(x).split(".")[0] for x in ids])


def main() -> None:
    h5ad_path = sys.argv[1]
    preds_path = sys.argv[2]
    out_path = sys.argv[3]
    overrides = parse_overrides(sys.argv[4:])

    adata = sc.read_h5ad(h5ad_path)

    if adata.raw is not None:
        counts = adata.raw.X
        if sp.issparse(counts):
            counts = counts.tocsr()
        else:
            counts = counts.copy()
        adata.raw = None
        adata.X = counts
    elif sp.issparse(adata.X):
        adata.X = adata.X.tocsr()

    adata.obsm.clear()
    adata.obsp.clear()
    adata.varm.clear()
    adata.varp.clear()
    adata.layers.clear()
    adata.uns.pop("log1p", None)
    adata.uns.pop("hvg", None)

    indiv_col = resolve_col(
        adata.obs, overrides.get("individual"), INDIV_ALIASES, "individual"
    )
    compartment_col = resolve_col(
        adata.obs, overrides.get("compartment"), COMPARTMENT_ALIASES, "compartment"
    )
    timepoint_col = resolve_col(
        adata.obs, overrides.get("timepoint"), TIMEPOINT_ALIASES, "timepoint"
    )
    batch_col = resolve_col(
        adata.obs, overrides.get("batch"), BATCH_ALIASES, "batch"
    )
    age_col = resolve_col(adata.obs, overrides.get("age"), AGE_ALIASES, "age")
    sex_col = resolve_col(adata.obs, overrides.get("sex"), SEX_ALIASES, "sex")

    adata.obs["individual"] = adata.obs[indiv_col].astype(str)
    adata.obs["compartment"] = normalise_compartment(adata.obs[compartment_col])
    adata.obs["timepoint"] = normalise_timepoint(adata.obs[timepoint_col])
    adata.obs["batch"] = adata.obs[batch_col].astype(str)
    adata.obs["age"] = pd.to_numeric(adata.obs[age_col], errors="coerce")
    adata.obs["sex"] = normalise_sex(adata.obs[sex_col])
    adata.obs["cell_label"] = (
        adata.obs["compartment"].astype(str) + " " + adata.obs["timepoint"].astype(str)
    )

    # Predictions CSV: barcode index / column + predicted_generation (+ optional pseudotime).
    preds = pd.read_csv(preds_path)
    barcode_col = None
    for cand in ["barcode", "cell_id", "cell", "obs_names", "Unnamed: 0"]:
        if cand in preds.columns:
            barcode_col = cand
            break
    if barcode_col is None:
        # First column is often the barcode index written without a name.
        barcode_col = preds.columns[0]
    preds = preds.rename(columns={barcode_col: "cell_id"})
    preds["cell_id"] = preds["cell_id"].astype(str)

    if "predicted_generation" not in preds.columns:
        raise SystemExit(
            f"predictions CSV missing predicted_generation; columns={list(preds.columns)}"
        )
    preds["division_state"] = preds["predicted_generation"]

    pt_pref = overrides.get("pseudotime")
    pt_col = None
    if pt_pref and pt_pref in preds.columns:
        pt_col = pt_pref
    elif pt_pref and pt_pref in adata.obs.columns:
        adata.obs["pseudotime"] = pd.to_numeric(adata.obs[pt_pref], errors="coerce")
    else:
        pt_col = first_col(preds, PSEUDOTIME_ALIASES)
        if pt_col is None:
            pt_col = first_col(adata.obs, PSEUDOTIME_ALIASES)
            if pt_col is not None:
                adata.obs["pseudotime"] = pd.to_numeric(
                    adata.obs[pt_col], errors="coerce"
                )
                pt_col = None

    keep_cols = ["cell_id", "division_state"]
    if pt_col is not None:
        preds = preds.rename(columns={pt_col: "pseudotime"})
        keep_cols.append("pseudotime")
    preds = preds[keep_cols].drop_duplicates(subset=["cell_id"])

    meta = adata.obs.copy()
    meta["cell_id"] = meta.index.astype(str)
    meta = meta.merge(preds, on="cell_id", how="left", suffixes=("", "_pred"))
    if "pseudotime_pred" in meta.columns:
        if "pseudotime" not in meta.columns:
            meta["pseudotime"] = meta["pseudotime_pred"]
        else:
            meta["pseudotime"] = meta["pseudotime"].fillna(meta["pseudotime_pred"])
        meta = meta.drop(columns=["pseudotime_pred"])
    if "pseudotime" not in meta.columns:
        raise SystemExit(
            "No pseudotime column found in h5ad or predictions CSV. "
            "Pass pseudotime=<col> after inspecting with bin/inspect-cytoimmgen.py"
        )

    meta = meta.set_index("cell_id")
    meta = meta.reindex(adata.obs_names.astype(str))
    adata.obs = meta
    adata.obs.index = adata.obs_names

    n_missing_div = int(adata.obs["division_state"].isna().sum())
    n_missing_pt = int(adata.obs["pseudotime"].isna().sum())
    if n_missing_div:
        print(f"Warning: {n_missing_div} cells missing division_state after join")
    if n_missing_pt:
        print(f"Warning: {n_missing_pt} cells missing pseudotime after join")

    symbol_col = first_col(adata.var, SYMBOL_ALIASES)
    if symbol_col is not None:
        adata.var["GeneSymbol"] = adata.var[symbol_col].astype(str)
    else:
        adata.var["GeneSymbol"] = adata.var_names.astype(str)

    adata.var_names = strip_ensembl_version(adata.var_names)
    adata.var_names_make_unique()

    # Drop cells that failed the predictions join for either interaction covariate.
    keep = adata.obs["division_state"].notna() & adata.obs["pseudotime"].notna()
    adata = adata[keep].copy()

    print(
        "Wrote harmonised h5ad with cell_label levels:",
        sorted(adata.obs["cell_label"].unique().tolist()),
    )
    adata.write_h5ad(out_path)


if __name__ == "__main__":
    main()
