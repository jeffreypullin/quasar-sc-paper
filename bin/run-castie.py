#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
import time
from pathlib import Path

SAMPLE_COVARS = "sex,age,PC_1,PC_2,geno_pc1,geno_pc2,geno_pc3,geno_pc4,geno_pc5,geno_pc6"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run CASTIE step1/step2/step3 per gene.")
    parser.add_argument("--cell-type", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--gene-list", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--subset-prefix", required=True)
    parser.add_argument("--chr-prefix", required=True)
    parser.add_argument("--anno", required=True)
    parser.add_argument(
        "--castie-dir",
        default="/home/jp2045/quasar-sc-paper/CASTIE",
    )
    parser.add_argument("--data-type", default="counts")
    parser.add_argument("--int-cov", default="pseudotime")
    parser.add_argument("--n-threads", type=int, default=8)
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
    return "-".join(str(text).strip().split())


def main() -> int:
    args = parse_args()
    genes = read_genes(args.gene_list)
    anno_by_gene = load_annotation(args.anno)
    env = build_env()
    pixi = [
        "pixi",
        "run",
        "--manifest-path",
        str(Path(args.castie_dir) / "pixi.toml"),
    ]

    gene_list_tag = safe_tag(Path(args.gene_list).name)
    data_type_tag = safe_tag(args.data_type)
    int_cov = args.int_cov
    all_covars = SAMPLE_COVARS + "," + int_cov
    manifest_path = Path(
        f"{safe_tag(args.cell_type)}-{safe_tag(args.chrom)}-{data_type_tag}-{gene_list_tag}-castie-files.tsv"
    )
    manifest_rows = []

    for gene in genes:
        output_prefix = f"{args.cell_type}-{gene}"
        step1_error = False
        variant_file = Path(f"{gene}_cis")

        start_time = time.time()
        step1_rc = run_cmd(
            pixi
            + [
                "Rscript",
                str(Path(args.castie_dir) / "extdata" / "step1_fitNULLGLMM_qtl.R"),
                "--useSparseGRMtoFitNULL=FALSE",
                "--useGRMtoFitNULL=FALSE",
                f"--phenoFile={args.input}",
                f"--phenoCol={gene}",
                f"--covarColList={all_covars}",
                f"--sampleCovarColList={SAMPLE_COVARS}",
                f"--dynamicCovarColList={int_cov}",
                "--sampleIDColinphenoFile=sample_id",
                "--traitType=count",
                f"--outputPrefix={output_prefix}-castie-step1",
                "--skipVarianceRatioEstimation=FALSE",
                "--isRemoveZerosinPheno=FALSE",
                "--isCovariateOffset=FALSE",
                "--isCovariateTransform=TRUE",
                "--skipModelFitting=FALSE",
                "--tol=0.00001",
                f"--plinkFile={args.subset_prefix}",
                "--IsOverwriteVarianceRatioFile=TRUE",
                "--offsetCol=log_cell_read_counts",
                "--tauInit=1,0.1,0",
                "--usePCG=FALSE",
                f"--nThreads={args.n_threads}",
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

        if not step1_error:
            chrom = args.chrom
            if isinstance(chrom, str) and chrom.lower().startswith("chr"):
                chrom = chrom[3:]
            write_region_file(anno_by_gene.get(gene, []))

            start_time = time.time()
            run_cmd(
                pixi
                + [
                    "Rscript",
                    str(Path(args.castie_dir) / "extdata" / "step2_tests_qtl.R"),
                    f"--bedFile={args.chr_prefix}.bed",
                    f"--bimFile={args.chr_prefix}.bim",
                    f"--famFile={args.chr_prefix}.fam",
                    f"--SAIGEOutputFile={gene}_cis",
                    f"--chrom={chrom}",
                    "--minMAF=0.05",
                    "--minMAC=5",
                    "--LOCO=FALSE",
                    f"--GMMATmodelFile={output_prefix}-castie-step1.rda",
                    "--SPAcutoff=2",
                    "--pval_cutoff_for_gxe=1",
                    f"--varianceRatioFile={output_prefix}-castie-step1.varianceRatio.txt",
                    "--rangestoIncludeFile=region-file.txt",
                    "--markers_per_chunk=1000",
                    "--output_format=txt",
                ],
                env,
            )
            end_time = time.time()
            step2_time = end_time - start_time

        manifest_rows.append(
            {
                "cell_type": args.cell_type,
                "chrom": args.chrom,
                "data_type": args.data_type,
                "int_cov": int_cov,
                "gene": gene,
                "variant_file": str(variant_file.resolve()) if variant_file.exists() else "NA",
                "step1_time": step1_time,
                "step2_time": step2_time,
            }
        )

    step3_time = 0
    region_file = Path("step3/step3_longformat.txt")
    if any(Path(row["variant_file"]).exists() for row in manifest_rows if row["variant_file"] != "NA"):
        start_time = time.time()
        run_cmd(
            pixi
            + [
                "python",
                str(Path(args.castie_dir) / "extdata" / "concat_step2_results.py"),
                "--input-dir=.",
                "--output=step3_input.txt",
                f"--contexts={int_cov}",
                "--file-pattern=*_cis",
                "--gene-regex=^(?P<gene>.+)_cis$",
                "--maf-min=0",
                "--maf-max=1",
            ],
            env,
        )
        run_cmd(
            pixi
            + [
                "Rscript",
                str(Path(args.castie_dir) / "extdata" / "step3_gene_pvalue.R"),
                "--input=step3_input.txt",
                "--outdir=step3",
            ],
            env,
        )
        step3_time = time.time() - start_time

    assigned_step3 = False
    for row in manifest_rows:
        if row["variant_file"] != "NA" and region_file.exists():
            row["region_file"] = str(region_file.resolve())
            row["step3_time"] = step3_time if not assigned_step3 else 0
            assigned_step3 = True
        else:
            row["region_file"] = "NA"
            row["step3_time"] = 0

    with open(manifest_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "cell_type",
                "chrom",
                "data_type",
                "int_cov",
                "gene",
                "variant_file",
                "region_file",
                "step1_time",
                "step2_time",
                "step3_time",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
