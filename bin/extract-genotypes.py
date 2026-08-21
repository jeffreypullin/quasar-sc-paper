#! /usr/bin/env python3
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

META_COLS = {"FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"}

plink_prefix = Path(sys.argv[1])
variants_tsv = Path(sys.argv[2])
out_path = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("genotype-dosages.tsv")

variants_df = pd.read_csv(variants_tsv, sep="\t")
if "snp_id" not in variants_df.columns:
    sys.exit(f"variants TSV must contain snp_id column; got: {', '.join(variants_df.columns)}")

VARIANTS = list(dict.fromkeys(variants_df["snp_id"].dropna().astype(str).tolist()))
if not VARIANTS:
    sys.exit("no snp_id values found in variants TSV")


def quasar_snp_id(chrom: str, pos: str, a1: str, a2: str) -> str:
    return f"{chrom}:{pos}{a2}-{a1}"


def norm_id(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", name).upper()


def find_variant_column(columns: list[str], plink_id: str, quasar_id: str) -> str | None:
    # plink2 --export A: columns are quasar_id + "_" + counted allele (e.g. 6:32608014T-C_T)
    export_cols = [c for c in columns if c.startswith(quasar_id + "_")]
    if len(export_cols) == 1:
        return export_cols[0]
    if quasar_id in columns:
        return quasar_id
    if plink_id in columns:
        return plink_id
    for col in columns:
        if col.startswith(plink_id + "_"):
            return col
    target = norm_id(plink_id)
    matches = [c for c in columns if norm_id(c) == target]
    if len(matches) == 1:
        return matches[0]
    return None


bim = pd.read_csv(
    f"{plink_prefix}.bim",
    sep="\t",
    header=None,
    names=["chrom", "snp", "cm", "pos", "a1", "a2"],
    dtype=str,
)
bim["quasar_id"] = bim.apply(
    lambda r: quasar_snp_id(r["chrom"], r["pos"], r["a1"], r["a2"]), axis=1
)

id_map = dict(zip(bim["quasar_id"], bim["snp"]))
missing = [v for v in VARIANTS if v not in id_map]
if missing:
    sys.exit(f"variants not in bim: {', '.join(missing)}")

extract_path = Path("extract-variants.snplist")
extract_path.write_text("\n".join(id_map[v] for v in VARIANTS) + "\n")

subprocess.run(
    [
        "plink2",
        "--bfile",
        str(plink_prefix),
        "--extract",
        str(extract_path),
        "--export",
        "A",
        "--out",
        "extracted",
        "--const-fid",
    ],
    check=True,
)

raw = pd.read_csv("extracted.raw", sep=r"\s+")
raw.columns = [c.lstrip("#") for c in raw.columns]
snp_cols = [c for c in raw.columns if c not in META_COLS]

rename = {}
for quasar_id in VARIANTS:
    plink_id = id_map[quasar_id]
    col = find_variant_column(snp_cols, plink_id, quasar_id)
    if col is None:
        sys.exit(
            f"no export column for {quasar_id} (plink id {plink_id}); "
            f"columns: {', '.join(snp_cols)}"
        )
    rename[col] = quasar_id

out = raw.rename(columns={"IID": "sample_id", **rename})
out = out[["sample_id"] + VARIANTS]
out.to_csv(out_path, sep="\t", index=False)
