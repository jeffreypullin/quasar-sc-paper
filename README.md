# Code accompanying the paper: "TODO"

Nextflow pipeline and analysis code for the paper describing the single-cell resolution and 
context-dependent eQTL mapping implemented in version 2.0 of [quasar](https://github.com/jeffreypullin/quasar).

## Run the pipeline

```bash
bash run.sh                          # default dataset (onek1k)
bash run.sh --datasets onek1k,tremor # choose datasets
bash run.sh --datasets randolph
```

## Datasets

Configured in `nextflow.config` under `params.dataset_configs` and selected with
`--datasets` (comma-separated; default `onek1k`):

| Dataset   | Description                                                                  |
|-----------|------------------------------------------------------------------------------|
| `onek1k`  | OneK1K PBMC cohort (14 cell types, GRCh37); main dataset for the paper.      |
| `tremor`  | Essential tremor cerebellum snRNA-seq (GRCh38); modules in `modules/tremor.nf`. |
| `randolph`| Influenza-infected PBMCs (Randolph et al.); expanded into `randolph_NI` and `randolph_flu` strata (`modules/randolph.nf`). |

