#! /usr/bin/env python3

import sys
from pathlib import Path

import scanpy as sc
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from cell_label_subset import cell_label_mask  # noqa: E402
import pooch
import zipfile

def flatten(xss):
    return [x for xs in xss for x in xs]

cell_label = sys.argv[1]
adata = sc.read_h5ad(sys.argv[2], backed="r")
subset = adata[
    cell_label_mask(adata.obs["cell_label"], cell_label).to_numpy()
].to_memory()

sc.pp.filter_genes(subset, min_cells=3)

subset.layers["counts"] = subset.X.copy()
sc.pp.normalize_total(subset)
sc.pp.log1p(subset)

# Mitochondrial percentaage.

subset.var["GeneSymbol"] = subset.var["GeneSymbol"].str.lstrip()

#mt_genes = sc.queries.mitochondrial_genes("hsapiens", attrname='ensembl_gene_id')
#subset.var["mt"] = subset.var_names.isin(mt_genes["ensembl_gene_id"])
subset.var["mt"] = subset.var["GeneSymbol"].str.startswith("MT-")
sc.pp.calculate_qc_metrics(subset, qc_vars=["mt"],  percent_top=None, inplace=True, log1p=True)

# Cell cycle.

p_zip = pooch.retrieve(
    "https://www.dropbox.com/s/3dby3bjsaf5arrw/cell_cycle_vignette_files.zip?dl=1",
    known_hash="sha256:6557fe2a2761d1550d41edf61b26c6d5ec9b42b4933d56c3161164a595f145ab",
    path="../data",
)

with zipfile.ZipFile(p_zip, "r") as f_zip:
    f_csv = zipfile.Path(f_zip, "nestorawa_forcellcycle_expressionMatrix.txt").open()
    cell_cycle_genes = zipfile.Path(f_zip, "regev_lab_cell_cycle_genes.txt").read_text().splitlines()

s_symbols = [x for x in cell_cycle_genes[:43] if x in subset.var["GeneSymbol"].values]
g2m_symbols = [x for x in cell_cycle_genes[43:] if x in subset.var["GeneSymbol"].values]
s_genes = pd.unique(subset.var_names[subset.var["GeneSymbol"].isin(s_symbols)]).tolist()
g2m_genes = pd.unique(subset.var_names[subset.var["GeneSymbol"].isin(g2m_symbols)]).tolist()

sc.tl.score_genes_cell_cycle(subset, s_genes=s_genes, g2m_genes=g2m_genes)

# Single-cell PCA.

sc.tl.pca(subset)

pc_data = subset.obsm['X_pca'][0:, 0:10]

pc_df = pd.DataFrame(
  data=pc_data[0:,0:],
  index=subset.obs_names,
  columns=[f"scPC_{i}" for i in range(1, 11)]
)

# Combine covariates.

subset_cov = pc_df.reset_index()
subset_cov["pct_counts_mt"] = subset.obs["pct_counts_mt"].values
subset_cov["S_score"] = subset.obs["S_score"].values
subset_cov["G2M_score"] = subset.obs["G2M_score"].values
subset_cov.insert(0, 'sample_id', subset.obs['individual'].values)
subset_cov.rename(columns={"barcode": "cell_id", "index": "cell_id"}, inplace=True)
subset_cov = subset_cov.sort_values(by='sample_id')

subset_cov.to_csv(f"{cell_label}-expr-covs.tsv", sep="\t", index=False)
