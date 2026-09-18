#!/usr/bin/env python3
# Emit the curated VCF IID → AnnData donor_id map, validated against a Tremor VCF.
#
# The published wells / Supp Table 1 order is NOT the GSA loading order. The map in
# resources/tremor-sample-map.tsv was recovered by matching published TensorQTL
# cis-eQTL dosages to pseudobulk expression (Hungarian assignment on INT-normalised
# slope-weighted fingerprints; LOCO stability ~86%; recovers TensorQTL effect signs).

import gzip
import sys
from pathlib import Path

import pandas as pd

map_path, vcf_path, out_path = sys.argv[1:4]

curated = pd.read_csv(map_path, sep="\t", dtype=str)
need = {"donor_id", "vcf_iid"}
missing = need - set(curated.columns)
if missing:
    raise SystemExit(f"curated map missing columns: {sorted(missing)}")

with gzip.open(vcf_path, "rt") as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            vcf_ids = line.rstrip("\n").split("\t")[9:]
            break
    else:
        raise SystemExit(f"no #CHROM line in {vcf_path}")

by_iid = curated.drop_duplicates("vcf_iid", keep="first").set_index("vcf_iid")
missing_ids = [iid for iid in vcf_ids if iid not in by_iid.index]
if missing_ids:
    raise SystemExit(
        f"{len(missing_ids)} VCF sample(s) absent from curated map "
        f"(e.g. {missing_ids[0]})"
    )

out = by_iid.loc[vcf_ids].reset_index()
# Keep donor_id first for plink --update-ids (newIID = column 1 in NF awk).
cols = ["donor_id", "vcf_iid"] + [
    c for c in out.columns if c not in {"donor_id", "vcf_iid"}
]
out[cols].to_csv(out_path, sep="\t", index=False)
