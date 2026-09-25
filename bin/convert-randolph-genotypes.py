#!/usr/bin/env python3
"""Convert Randolph MatrixEQTL genotype dosages to per-chromosome PLINK beds.

The supplied genotypes.txt is the authors' QC-filtered set (MAF > 5%,
missingness < 10%, HWE P >= 1e-5; 6,305,923 SNPs). Missing dosages = -9.

Usage:
  convert-randolph-genotypes.py <genotypes.txt> <SNP_positions.txt> <out_prefix>

Writes <out_prefix>-chr{1..22}.{bed,bim,fam}. Requires plink2 on PATH.
"""

import gzip
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


def parse_alleles(snp):
    # IDs look like 1_13110_G_A
    _chrom, _pos, a1, a2 = snp.split("_", 3)
    return a1, a2


def main():
    if len(sys.argv) < 4:
        raise SystemExit(
            "Usage: convert-randolph-genotypes.py "
            "<genotypes.txt> <SNP_positions.txt> <out_prefix>"
        )
    geno_path = Path(sys.argv[1])
    snp_path = Path(sys.argv[2])
    out_prefix = sys.argv[3]

    snp_chr = {}
    snp_pos = {}
    with open(snp_path) as fh:
        next(fh)
        for line in fh:
            parts = line.replace('"', "").split()
            if len(parts) < 3:
                continue
            # Observed: idx snp chr pos
            if len(parts) >= 4 and parts[0].isdigit():
                snp, chrom, pos = parts[1], parts[2], parts[3]
            else:
                snp, chrom, pos = parts[0], parts[1], parts[2]
            chrom = chrom.replace("chr", "")
            if chrom not in {str(i) for i in range(1, 23)}:
                continue
            snp_chr[snp] = chrom
            snp_pos[snp] = int(pos)

    tmpdir = Path(tempfile.mkdtemp(prefix="randolph-geno-"))
    writers = {}
    counts = defaultdict(int)
    samples = []

    def writer_for(chrom):
        if chrom not in writers:
            path = tmpdir / f"chr{chrom}.dose.gz"
            # noheader dosage: SNP A1 A2 <doses in .fam order>
            writers[chrom] = gzip.open(path, "wt")
        return writers[chrom]

    with open(geno_path) as fh:
        # MatrixEQTL-style: header is sample IDs only (no SNP column name).
        header = fh.readline().rstrip("\n").split("\t")
        samples = [s.strip('"') for s in header]
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            snp = parts[0].strip('"')
            if snp not in snp_chr:
                continue
            chrom = snp_chr[snp]
            try:
                a1, a2 = parse_alleles(snp)
            except ValueError:
                continue
            dosages = [
                "NA" if x.strip() in ("-9", "NA", "nan", "") else x.strip()
                for x in parts[1:]
            ]
            if len(dosages) != len(samples):
                raise SystemExit(
                    f"{snp}: expected {len(samples)} dosages, got {len(dosages)}"
                )
            w = writer_for(chrom)
            w.write(f"{snp}\t{a1}\t{a2}\t")
            w.write("\t".join(dosages))
            w.write("\n")
            counts[chrom] += 1

    for w in writers.values():
        w.close()

    for chrom in counts:
        map_path = tmpdir / f"chr{chrom}.map"
        with gzip.open(tmpdir / f"chr{chrom}.dose.gz", "rt") as din, open(
            map_path, "w"
        ) as mout:
            for line in din:
                snp = line.split("\t", 1)[0]
                mout.write(f"{chrom} {snp} 0 {snp_pos[snp]}\n")

    fam = tmpdir / "samples.fam"
    with open(fam, "w") as sf:
        for s in samples:
            # FID IID father mother sex phenotype
            sf.write(f"0\t{s}\t0\t0\t1\t-9\n")

    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    for chrom in sorted(counts, key=int):
        out = f"{out_prefix}-chr{chrom}"
        # Header sample tokens are IIDs only; plink2 would parse them as
        # FID/IID pairs. Use noheader + .fam order instead.
        cmd = [
            "plink2",
            "--import-dosage",
            str(tmpdir / f"chr{chrom}.dose.gz"),
            "noheader",
            "format=1",
            "ref-first",
            "--map",
            str(tmpdir / f"chr{chrom}.map"),
            "--fam",
            str(fam),
            "--make-bed",
            "--out",
            out,
            "--output-chr",
            "26",
            "--set-all-var-ids",
            "@:#$r-$a",
        ]
        print("Running:", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
        print(f"chr{chrom}: {counts[chrom]} variants → {out}.bed", flush=True)

    for p in tmpdir.glob("*"):
        try:
            p.unlink()
        except OSError:
            pass
    try:
        tmpdir.rmdir()
    except OSError:
        pass


if __name__ == "__main__":
    main()
