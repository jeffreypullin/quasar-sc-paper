#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
import time
from pathlib import Path

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run SAIGE-QTL step1/step2 per gene.")
    parser.add_argument("--cell-type", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--gene-list", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--subset-prefix", required=True)
    parser.add_argument("--chr-prefix", required=True)
    parser.add_argument("--anno", required=True)
    return parser.parse_args()

def build_env() -> dict:
    env = os.environ.copy()
    env["R_LIBS_USER"] = ""
    env["R_PROFILE_USER"] = "/dev/null"
    env["R_ENVIRON_USER"] = "/dev/null"
    return env

def run_cmd(args: list, env: dict) -> int:
    return subprocess.run(args, env=env).returncode

def read_genes(path: str) -> list:
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]

def load_annotation(path: str) -> dict:
    by_gene = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row.get("phenotype_id")
            if not gene:
                continue
            by_gene.setdefault(gene, []).append(row)
    return by_gene

def write_region_file(entries: list) -> None:
    lines = []
    for row in entries:
        chrom = row.get("#chr")
        start = int(float(row["start"]))
        end = int(float(row["end"]))
        lines.append(f"{chrom}\t{start - 1000000}\t{end + 1000000}")
    with open("region-file.txt", "w", encoding="utf-8") as file:
        file.write("\n".join(lines))

def safe_tag(text: str) -> str:
    # Keep filenames simple and stable.
    return "-".join(str(text).strip().split())

def main() -> int:
    args = parse_args()
    genes = read_genes(args.gene_list)
    anno_by_gene = load_annotation(args.anno)
    env = build_env()

    gene_list_tag = safe_tag(Path(args.gene_list).name)
    manifest_path = Path(f"{safe_tag(args.cell_type)}-{safe_tag(args.chrom)}-{gene_list_tag}-saigeqtl-files.tsv")
    manifest_rows = []

    for gene in genes:

        output_prefix = f"{args.cell_type}-{gene}"
        step1_error = False
        variant_file = Path(f"{output_prefix}-variant")
        region_file = Path(f"{output_prefix}-region")
        
        start_time = time.time()
        step1_rc = run_cmd(
            [
                "pixi",
                "run",
                "--manifest-path",
                "/home/jp2045/quasar-sc-paper/SAIGEQTL/pixi.toml",
                "Rscript",
                "/home/jp2045/quasar-sc-paper/SAIGEQTL/extdata/step1_fitNULLGLMM_qtl.R",
                "--useSparseGRMtoFitNULL=FALSE",
                "--useGRMtoFitNULL=FALSE",
                f"--phenoFile={args.input}",
                f"--phenoCol={gene}",
                "--covarColList=sex,age,PC_1,PC_2,geno_pc1,geno_pc2,geno_pc3,geno_pc4,geno_pc5,geno_pc6",
                "--sampleCovarColList=sex,age,PC_1,PC_2,geno_pc1,geno_pc2,geno_pc3,geno_pc4,geno_pc5,geno_pc6",
                "--sampleIDColinphenoFile=sample_id",
                "--traitType=count",
                f"--outputPrefix={output_prefix}-saigeqtl-step1",
                "--skipVarianceRatioEstimation=FALSE",
                "--isRemoveZerosinPheno=FALSE",
                "--isCovariateOffset=FALSE",
                "--isCovariateTransform=TRUE",
                "--skipModelFitting=FALSE",
                "--tol=0.00001",
                f"--plinkFile={args.subset_prefix}",
                "--IsOverwriteVarianceRatioFile=TRUE",
            ],
            env,
        )
        end_time = time.time()

        if step1_rc == 0:
            step1_time = end_time - start_time
        else:
            step1_error = True
            step1_time = 0
            step2_time = 0
            step3_time = 0

        if not step1_error:

            # Convert chromosome labels to number.
            chrom = args.chrom
            if isinstance(chrom, str) and chrom.lower().startswith("chr"):
                chrom = chrom[3:]
            write_region_file(anno_by_gene.get(gene, []))
            
            start_time = time.time()
            run_cmd(
                [
                    "pixi",
                    "run",
                    "--manifest-path",
                    "/home/jp2045/quasar-sc-paper/SAIGEQTL/pixi.toml",
                    "Rscript",
                    "/home/jp2045/quasar-sc-paper/SAIGEQTL/extdata/step2_tests_qtl.R",
                    f"--bedFile={args.chr_prefix}.bed",
                    f"--bimFile={args.chr_prefix}.bim",
                    f"--famFile={args.chr_prefix}.fam",
                    f"--SAIGEOutputFile={output_prefix}-variant",
                    f"--chrom={chrom}",
                    "--minMAF=0",
                    "--minMAC=1",
                    "--LOCO=FALSE",
                    f"--GMMATmodelFile={output_prefix}-saigeqtl-step1.rda",
                    "--SPAcutoff=2",
                    f"--varianceRatioFile={output_prefix}-saigeqtl-step1.varianceRatio.txt",
                    "--rangestoIncludeFile=region-file.txt",
                    "--markers_per_chunk=10000",
                ],
                env,
            )
            end_time = time.time()
            step2_time = end_time - start_time

        if not step1_error:
            start_time = time.time()
            run_cmd(
                [
                "pixi",
                "run",
                "--manifest-path",
                "/home/jp2045/quasar-sc-paper/SAIGEQTL/pixi.toml",
                "Rscript",
                "/home/jp2045/quasar-sc-paper/SAIGEQTL/extdata/step3_gene_pvalue_qtl.R",
                f"--assocFile={output_prefix}-variant",
                f"--geneName={gene}",
                f"--genePval_outputFile={output_prefix}-region",
                ],
                env,
            )
            end_time = time.time()
            step3_time = end_time - start_time

        manifest_rows.append(
            {
                "cell_type": args.cell_type,
                "chrom": args.chrom,
                "gene": gene,
                "variant_file": str(variant_file.resolve()) if variant_file.exists() else "NA",
                "region_file": str(region_file.resolve()) if region_file.exists() else "NA",
                "step1_time": step1_time,
                "step2_time": step2_time,
                "step3_time": step3_time
            }
        )

    with open(manifest_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "cell_type",
                "chrom",
                "gene",
                "variant_file",
                "region_file",
                "step1_time",
                "step2_time",
                "step3_time"
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
