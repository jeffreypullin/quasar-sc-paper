#!/usr/bin/env python3
"""Build QuASAR annotation BED from Randolph GRCh38.92 gene positions.

Usage:
  convert-randolph-annot.py <GRCh38.92_gene_positions.txt> <out.bed>

Output columns: #chr start end phenotype_id
Uses gene-body start (S1) as a point-like TSS for QuASAR windowing,
matching the gene-symbol IDs used as Seurat features.
"""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(
            "Usage: convert-randolph-annot.py <gene_positions.txt> <out.bed>"
        )
    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    rows: list[tuple[int, int, int, str]] = []
    with open(in_path) as fh:
        next(fh)
        for line in fh:
            # File is inconsistently quoted; strip all quotes then split.
            parts = line.replace('"', "").split()
            # Observed: idx Gene_ID chromosome S1 S2
            if len(parts) < 5:
                continue
            gene, chrom, start, end = parts[1], parts[2], parts[3], parts[4]
            chrom = chrom.replace("chr", "")
            if chrom not in {str(i) for i in range(1, 23)}:
                continue
            start_i, end_i = int(start), int(end)
            if end_i < start_i:
                start_i, end_i = end_i, start_i
            tss = start_i
            rows.append((int(chrom), tss, tss + 1, gene))

    rows.sort()
    with open(out_path, "w") as out:
        out.write("#chr\tstart\tend\tphenotype_id\n")
        for chrom, start, end, gene in rows:
            out.write(f"{chrom}\t{start}\t{end}\t{gene}\n")
    print(f"Wrote {len(rows)} autosomal genes to {out_path}")


if __name__ == "__main__":
    main()
