#!/usr/bin/env python3

import sys
from pathlib import Path

import pandas as pd
import scanpy as sc
import scipy.sparse as sp

CELL_LABEL_COLS = ["author_cell_type", "cell_label", "cell_type", "annotation"]
INDIV_COLS = ["individual", "donor_id", "donor", "sample_id"]
SEX_COLS = ["sex", "Sex"]
AGE_COLS = ["age", "Age"]
SYMBOL_COLS = ["GeneSymbol", "feature_name", "gene_symbols", "name"]

# Match OneK1K / PLINK coding: 1 = male, 2 = female.
SEX_ONTOLOGY_MAP = {
    "PATO:0000384": 1,  # male
    "PATO:0000383": 2,  # female
}
DONOR_TYPOS = {"SK00397": "SK00387"}
CELL_LABEL_MAP = {
    "granule cell": "Granule",
    "granule cells": "Granule",
    "Granule cells": "Granule",
    "Purkinje cell": "Purkinje",
    "Purkinje cells": "Purkinje",
    "molecular layer interneuron": "MLI 1",
    "Bergmann glial cell": "Bergmann",
    "Bergmann glia": "Bergmann",
    "astrocyte": "Astrocytes",
    "microglial cell": "Microglia",
    "endothelial cell": "Endocytes",
    "endocyte": "Endocytes",
    "endocytes": "Endocytes",
    "pericyte": "Pericytes",
    "oligodendrocyte": "Oligodendrocytes",
    "oligodendrocyte precursor cell": "OPC",
    "Golgi cell": "Golgi",
    "Golgi cells": "Golgi",
    "unipolar brush cell": "UBC",
    "unipolar brush cells": "UBC",
}


def first_col(frame, names):
    for name in names:
        if name in frame.columns:
            return name
    return None


adata = sc.read_h5ad(sys.argv[1])
out_path = sys.argv[2] if len(sys.argv) > 2 else "harmonised.h5ad"
sample_map_path = sys.argv[3] if len(sys.argv) > 3 else None

# Tremor ships log-normalised values in X and integer counts in raw.X.
# Keep counts only so downstream scripts see the same matrix shape as OneK1K.
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

label_col = first_col(adata.obs, CELL_LABEL_COLS)
adata.obs["cell_label"] = (
    adata.obs[label_col]
    .astype(str)
    .replace(CELL_LABEL_MAP)
    .str.replace("_", " ", regex=False)
)

indiv_col = first_col(adata.obs, INDIV_COLS)
adata.obs["individual"] = (
    adata.obs[indiv_col].astype(str).replace(DONOR_TYPOS)
)

if sample_map_path is not None:
    keep = (
        pd.read_csv(sample_map_path, sep="\t")["donor_id"].astype(str).unique().tolist()
    )
    adata = adata[adata.obs["individual"].isin(keep)].copy()

sex_col = first_col(adata.obs, SEX_COLS)
if sex_col is not None:
    adata.obs["sex"] = adata.obs[sex_col]
elif "sex_ontology_term_id" in adata.obs.columns:
    adata.obs["sex"] = (
        adata.obs["sex_ontology_term_id"].astype(str).map(SEX_ONTOLOGY_MAP).astype(int)
    )

age_col = first_col(adata.obs, AGE_COLS)
if age_col is not None:
    adata.obs["age"] = adata.obs[age_col]

symbol_col = first_col(adata.var, SYMBOL_COLS)
if symbol_col is not None:
    adata.var["GeneSymbol"] = adata.var[symbol_col].astype(str)
else:
    adata.var["GeneSymbol"] = adata.var_names.astype(str)

adata.write_h5ad(out_path)
