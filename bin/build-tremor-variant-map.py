#!/usr/bin/env python3
"""Map original VCF variant IDs to QuASAR/plink @:#$r-$a IDs."""

import gzip
import sys
from pathlib import Path

out_path = sys.argv[1]
vcf_paths = sys.argv[2:]


def norm_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.startswith("chr") else chrom


with open(out_path, "w") as out:
    out.write("variant_id\tchrom\tpos\tref\talt\tquasar_snp_id\n")
    for vcf_path in vcf_paths:
        opener = gzip.open if str(vcf_path).endswith(".gz") else open
        with opener(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                chrom, pos, vid, ref, alt, *_ = line.rstrip("\n").split("\t")
                if alt == "." or "," in alt:
                    continue
                chrom_n = norm_chrom(chrom)
                quasar_id = f"{chrom_n}:{pos}{ref}-{alt}"
                out.write(
                    f"{vid}\t{chrom_n}\t{pos}\t{ref}\t{alt}\t{quasar_id}\n"
                )
