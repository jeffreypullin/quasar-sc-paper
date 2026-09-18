#!/usr/bin/env python3
"""Harmonise Tremor QuASAR LM and published TensorQTL nominal results."""

import argparse
import hashlib
import zipfile

import numpy as np
import pandas as pd


def norm_chrom(chrom: str) -> str:
    chrom = str(chrom)
    return chrom[3:] if chrom.startswith("chr") else chrom


def pearson(x, y):
    if len(x) < 2:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x, y):
    if len(x) < 2:
        return np.nan
    xr = pd.Series(x).rank().to_numpy()
    yr = pd.Series(y).rank().to_numpy()
    return pearson(xr, yr)


def load_gene_map(path):
    gene_map = pd.read_csv(path, sep="\t", dtype=str)
    if "gene_symbol" not in gene_map.columns:
        for alt in ("GeneSymbol", "name"):
            if alt in gene_map.columns:
                gene_map = gene_map.rename(columns={alt: "gene_symbol"})
                break
    return gene_map[["feature_id", "gene_symbol"]].astype(str)


def load_annot(path):
    annot = pd.read_csv(path, sep="\t", dtype={"phenotype_id": str})
    annot["chrom"] = annot["#chr"].map(norm_chrom)
    annot["tss"] = annot["start"].astype(int)
    return annot[["phenotype_id", "chrom", "tss"]].rename(
        columns={"phenotype_id": "feature_id"}
    )


def resolve_symbols(gene_map, annot):
    """Keep unambiguous Ensembl↔symbol links; drop multi-mapped symbols."""
    gm = gene_map.merge(annot, on="feature_id", how="left")
    symbol_n = gm.groupby("gene_symbol")["feature_id"].nunique()
    ambiguous_symbols = set(symbol_n[symbol_n > 1].index)
    unique = gm[~gm["gene_symbol"].isin(ambiguous_symbols)].copy()
    return unique, len(ambiguous_symbols)


def load_variant_map(path):
    vmap = pd.read_csv(
        path,
        sep="\t",
        dtype={
            "variant_id": str,
            "chrom": str,
            "ref": str,
            "alt": str,
            "quasar_snp_id": str,
        },
    )
    vmap["chrom"] = vmap["chrom"].map(norm_chrom)
    vmap["pos"] = vmap["pos"].astype(int)
    return vmap.drop_duplicates("variant_id", keep="first")


def load_quasar_variants(paths, gene_resolved):
    pieces = []
    for path in paths:
        df = pd.read_csv(
            path,
            sep="\t",
            usecols=[
                "feature_id",
                "snp_id",
                "chrom",
                "pos",
                "alt",
                "ref",
                "beta",
                "se",
                "pvalue",
            ],
            dtype={
                "feature_id": str,
                "snp_id": str,
                "chrom": str,
                "alt": str,
                "ref": str,
            },
        )
        df["chrom"] = df["chrom"].map(norm_chrom)
        df["pos"] = df["pos"].astype(int)
        pieces.append(df)
    quasar = pd.concat(pieces, ignore_index=True)
    quasar = quasar.merge(
        gene_resolved[["feature_id", "gene_symbol"]],
        on="feature_id",
        how="inner",
    )
    quasar = quasar.rename(
        columns={
            "beta": "beta_quasar",
            "se": "se_quasar",
            "pvalue": "pvalue_quasar",
            "snp_id": "quasar_snp_id",
        }
    )
    quasar["z_quasar"] = quasar["beta_quasar"] / quasar["se_quasar"]
    return quasar


def iter_tensorqtl_parquets(zip_path):
    with zipfile.ZipFile(zip_path) as zf:
        names = sorted(
            n
            for n in zf.namelist()
            if n.endswith(".parquet") and not n.startswith("__")
        )
        for name in names:
            with zf.open(name) as fh:
                yield name, pd.read_parquet(fh)


def sample_rows(df, n, seed_key):
    """Uniform random subsample for plotting."""
    if len(df) <= n:
        return df.copy()
    seed = int(hashlib.md5(seed_key.encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(df), size=n, replace=False)
    return df.iloc[np.sort(idx)].copy()


def stratified_z_metrics(z_q, z_t, abs_thresh):
    mask = np.abs(z_t) >= abs_thresh
    if mask.sum() < 2:
        return np.nan, np.nan, 0
    return (
        pearson(z_q[mask], z_t[mask]),
        float(np.mean(np.sign(z_q[mask]) == np.sign(z_t[mask]))),
        int(mask.sum()),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-type", required=True)
    parser.add_argument("--tensorqtl-zip", required=True)
    parser.add_argument("--gene-map", required=True)
    parser.add_argument("--variant-map", required=True)
    parser.add_argument("--annot-bed", required=True)
    parser.add_argument("--quasar-variants", nargs="+", required=True)
    parser.add_argument("--sample-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--sample-n", type=int, default=10000)
    args = parser.parse_args()

    gene_map = load_gene_map(args.gene_map)
    annot = load_annot(args.annot_bed)
    gene_resolved, n_ambiguous_symbols = resolve_symbols(gene_map, annot)
    vmap = load_variant_map(args.variant_map)
    vmap_by_id = vmap.set_index("variant_id", drop=False)

    quasar = load_quasar_variants(args.quasar_variants, gene_resolved)
    n_quasar = len(quasar)
    q_symbols = set(quasar["gene_symbol"].unique())

    matched_chunks = []
    n_tq = 0
    n_tq_mapped = 0
    n_allele_mismatch = 0
    n_unmatched_variant = 0
    n_unmatched_gene = 0

    for _name, tq in iter_tensorqtl_parquets(args.tensorqtl_zip):
        n_tq += len(tq)
        keep_cols = [
            c
            for c in (
                "phenotype_id",
                "variant_id",
                "pval_nominal",
                "slope",
                "slope_se",
            )
            if c in tq.columns
        ]
        tq = tq[keep_cols].copy()
        tq["phenotype_id"] = tq["phenotype_id"].astype(str)
        tq["variant_id"] = tq["variant_id"].astype(str)

        mapped = tq.join(vmap_by_id, on="variant_id", how="left", rsuffix="_map")
        missing_var = mapped["quasar_snp_id"].isna()
        n_unmatched_variant += int(missing_var.sum())
        mapped = mapped.loc[~missing_var].copy()
        n_tq_mapped += len(mapped)

        gene_ok = mapped["phenotype_id"].isin(q_symbols)
        n_unmatched_gene += int((~gene_ok).sum())
        mapped = mapped.loc[gene_ok].copy()
        if mapped.empty:
            continue

        mapped = mapped.rename(
            columns={
                "phenotype_id": "gene_symbol",
                "pval_nominal": "pvalue_tensorqtl",
                "slope": "beta_tensorqtl",
                "slope_se": "se_tensorqtl",
                "ref": "ref_tq",
                "alt": "alt_tq",
            }
        )
        mapped["chrom"] = mapped["chrom"].map(norm_chrom)
        mapped["pos"] = mapped["pos"].astype(int)
        mapped = mapped[
            [
                "gene_symbol",
                "variant_id",
                "chrom",
                "pos",
                "ref_tq",
                "alt_tq",
                "pvalue_tensorqtl",
                "beta_tensorqtl",
                "se_tensorqtl",
            ]
        ]

        pos_join = mapped.merge(
            quasar,
            on=["gene_symbol", "chrom", "pos"],
            how="inner",
        )
        if pos_join.empty:
            continue

        flip = np.where(
            (pos_join["ref_tq"] == pos_join["ref"])
            & (pos_join["alt_tq"] == pos_join["alt"]),
            1,
            np.where(
                (pos_join["ref_tq"] == pos_join["alt"])
                & (pos_join["alt_tq"] == pos_join["ref"]),
                -1,
                0,
            ),
        )
        n_allele_mismatch += int((flip == 0).sum())
        keep = flip != 0
        pos_join = pos_join.loc[keep].copy()
        pos_join["allele_flip"] = flip[keep]
        pos_join["beta_tensorqtl"] = (
            pos_join["beta_tensorqtl"] * pos_join["allele_flip"]
        )
        pos_join["z_tensorqtl"] = (
            pos_join["beta_tensorqtl"] / pos_join["se_tensorqtl"]
        )
        matched_chunks.append(
            pos_join[
                [
                    "gene_symbol",
                    "feature_id",
                    "variant_id",
                    "quasar_snp_id",
                    "chrom",
                    "pos",
                    "ref",
                    "alt",
                    "allele_flip",
                    "beta_quasar",
                    "se_quasar",
                    "pvalue_quasar",
                    "z_quasar",
                    "beta_tensorqtl",
                    "se_tensorqtl",
                    "pvalue_tensorqtl",
                    "z_tensorqtl",
                ]
            ]
        )

    if matched_chunks:
        matched = pd.concat(matched_chunks, ignore_index=True)
        matched = matched.drop_duplicates(
            ["feature_id", "quasar_snp_id"], keep="first"
        )
    else:
        matched = pd.DataFrame(
            columns=[
                "gene_symbol",
                "feature_id",
                "variant_id",
                "quasar_snp_id",
                "chrom",
                "pos",
                "ref",
                "alt",
                "allele_flip",
                "beta_quasar",
                "se_quasar",
                "pvalue_quasar",
                "z_quasar",
                "beta_tensorqtl",
                "se_tensorqtl",
                "pvalue_tensorqtl",
                "z_tensorqtl",
            ]
        )

    finite = matched[
        np.isfinite(matched["z_quasar"])
        & np.isfinite(matched["z_tensorqtl"])
        & np.isfinite(matched["pvalue_quasar"])
        & np.isfinite(matched["pvalue_tensorqtl"])
        & (matched["se_quasar"] > 0)
        & (matched["se_tensorqtl"] > 0)
    ].copy()

    sign_agree = (
        float(np.mean(np.sign(finite["z_quasar"]) == np.sign(finite["z_tensorqtl"])))
        if len(finite)
        else np.nan
    )
    z_q = finite["z_quasar"].to_numpy()
    z_t = finite["z_tensorqtl"].to_numpy()
    p_q = -np.log10(np.clip(finite["pvalue_quasar"].to_numpy(), 1e-300, 1))
    p_t = -np.log10(np.clip(finite["pvalue_tensorqtl"].to_numpy(), 1e-300, 1))
    pearson_z2, sign_z2, n_z2 = stratified_z_metrics(z_q, z_t, 2.0)
    pearson_z3, sign_z3, n_z3 = stratified_z_metrics(z_q, z_t, 3.0)

    summary = pd.DataFrame(
        [
            {
                "cell_type": args.cell_type,
                "n_quasar_tests": n_quasar,
                "n_tensorqtl_tests": n_tq,
                "n_tensorqtl_variant_mapped": n_tq_mapped,
                "n_tensorqtl_unmatched_variant": n_unmatched_variant,
                "n_tensorqtl_unmatched_gene": n_unmatched_gene,
                "n_ambiguous_symbols_dropped": n_ambiguous_symbols,
                "n_allele_mismatch": n_allele_mismatch,
                "n_matched": len(matched),
                "n_finite": len(finite),
                "pearson_z": pearson(z_q, z_t),
                "spearman_z": spearman(z_q, z_t),
                "pearson_neglog10p": pearson(p_q, p_t),
                "spearman_neglog10p": spearman(p_q, p_t),
                "sign_concordance": sign_agree,
                "n_flipped_alleles": int((finite["allele_flip"] == -1).sum())
                if len(finite)
                else 0,
            }
        ]
    )

    sample = sample_rows(finite, args.sample_n, args.cell_type)
    sample["cell_type"] = args.cell_type
    sample["neglog10p_quasar"] = -np.log10(
        np.clip(sample["pvalue_quasar"], 1e-300, 1)
    )
    sample["neglog10p_tensorqtl"] = -np.log10(
        np.clip(sample["pvalue_tensorqtl"], 1e-300, 1)
    )

    sample.to_csv(args.sample_out, sep="\t", index=False)
    summary.to_csv(args.summary_out, sep="\t", index=False)


if __name__ == "__main__":
    main()
