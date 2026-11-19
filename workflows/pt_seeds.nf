include { RUN_SLINGSHOT_SEED ; JOIN_SC_INT_COVS as JOIN_SC_INT_COVS_SEED } from '../modules/preprocess'
include { COMBINE_SC_COVS as COMBINE_SC_COVS_SEED ; FILTER_COVS as FILTER_COVS_SEED ;
          RUN_QUASAR_SC as RUN_QUASAR_SC_PT_SEED } from '../modules/quasar'
include { COMPUTE_CELL_GROUPS as COMPUTE_CELL_GROUPS_PT_SEED } from '../modules/cellgroups'
include { PLOT_SLINGSHOT_SEED ; TABULATE_PT_SEED_EGENES } from '../modules/analysis'

workflow PT_SEEDS {
    take:
    slingshot_adata
    sc_expr_covs
    pb_expr_covs
    geno_pcs
    sc_counts
    anno_bed
    bed_files
    sc_covs
    gene_map

    main:
    // Extra PHATE seeds. The published trajectory is the default seed inside
    // run-slingshot.R (20261013) and is not passed through here.
    pt_seeds = Channel.of(1, 2, 3, 4, 5)
    slingshot_seed_out = RUN_SLINGSHOT_SEED(slingshot_adata.combine(pt_seeds))
    PLOT_SLINGSHOT_SEED(
        slingshot_seed_out.map { info, seed, _pt, plot_data -> tuple(info, seed, plot_data) }
    )

    // Same OneK1K B_all counts model as the published run: P-GLMM GxC
    // (k = none, -i pseudotime) and K = 5 grouped bins. int_cov is seed-tagged
    // so these tasks stay out of the published pseudotime figures. Grouped
    // fits still use the donor-level bulk_pca covariates; pseudotime only
    // defines the bins.
    def pt_seed_int_cov = { seed -> "pseudotime_seed${seed}".toString() }
    pt_seed_expr = sc_expr_covs
        .filter { it[0].dataset == "onek1k" && it[0].cell_type == "B_all" && it[0].is_full_data }
        .map { info, expr -> tuple(info.dataset, info, expr) }
        .combine(
            slingshot_seed_out.map { sinfo, seed, pt, _plot -> tuple(sinfo.dataset, seed, pt) },
            by: 0
        )
        .map { dataset, info, expr, seed, pt -> tuple(info + [pt_seed: seed], expr, pt) }
    pt_seed_expr_aug = JOIN_SC_INT_COVS_SEED(pt_seed_expr)
    pt_seed_pb_expr = pb_expr_covs
        .filter { it[0].dataset == "onek1k" && it[0].cell_type == "B_all" && it[0].is_full_data }
        .map { info, pb -> tuple(info.dataset, pb) }
    pt_seed_covs = FILTER_COVS_SEED(
        COMBINE_SC_COVS_SEED(
            pt_seed_expr_aug
                .map { info, sc_c -> tuple(info.dataset, info, sc_c) }
                .combine(pt_seed_pb_expr, by: 0)
                .combine(geno_pcs, by: 0)
                .map { dataset, info, sc_c, pb_c, pcs -> tuple(info, pb_c, sc_c, pcs) }
        ).map { info, cov -> tuple(info, cov, "bulk_pca+pseudotime") }
    )
    pt_seed_pheno = sc_counts
        .filter { _dt, info, _pheno ->
            info.dataset == "onek1k" && info.cell_type == "B_all" && info.is_full_data
        }
        .map { _dt, info, pheno -> tuple(info.dataset, pheno) }
    pt_seed_anno = anno_bed
        .filter { dataset, _bed -> dataset == "onek1k" }
    pt_seed_beds = bed_files
        .filter { dataset, _chr, _bed -> dataset == "onek1k" }
    pt_seed_bulk_pca = sc_covs
        .filter { info, cov_spec, _covs ->
            cov_spec == "bulk_pca" &&
                info.dataset == "onek1k" &&
                info.cell_type == "B_all" &&
                info.is_full_data
        }
        .map { info, _cov_spec, covs -> tuple(info.dataset, covs) }
    pt_seed_gxc = pt_seed_covs
        .map { info, _cov_spec, covs -> tuple(info.dataset, info, covs) }
        .combine(pt_seed_pheno, by: 0)
        .combine(pt_seed_anno, by: 0)
        .combine(pt_seed_beds, by: 0)
        .map { dataset, info, covs, pheno, anno, chr, bed ->
            tuple(
                info + [
                    chr: chr,
                    cov_spec: "bulk_pca+pseudotime",
                    data_type: "counts",
                    model: "p_glmm_sc",
                    int_cov: pt_seed_int_cov(info.pt_seed),
                    interaction_cov: "pseudotime",
                    k: "none",
                ],
                pheno, covs, anno, bed, covs
            )
        }
    pt_seed_groups = COMPUTE_CELL_GROUPS_PT_SEED(
        pt_seed_covs.map { info, _cov_spec, covs ->
            tuple(
                info + [
                    int_cov: pt_seed_int_cov(info.pt_seed),
                    group_cov: "pseudotime",
                    k: 5,
                ],
                covs
            )
        }
    )
    pt_seed_grouped = pt_seed_groups
        .map { info, groups -> tuple(info.dataset, info, groups) }
        .combine(pt_seed_pheno, by: 0)
        .combine(pt_seed_bulk_pca, by: 0)
        .combine(pt_seed_anno, by: 0)
        .combine(pt_seed_beds, by: 0)
        .map { dataset, info, groups, pheno, bulk_covs, anno, chr, bed ->
            tuple(
                info + [
                    chr: chr,
                    cov_spec: "bulk_pca+pseudotime",
                    data_type: "counts",
                    model: "p_glmm_sc",
                    int_cov: pt_seed_int_cov(info.pt_seed),
                    interaction_cov: "pseudotime",
                    k: 5,
                ],
                pheno, bulk_covs, anno, bed, groups
            )
        }
    pt_seed_quasar = RUN_QUASAR_SC_PT_SEED(pt_seed_gxc.mix(pt_seed_grouped))
    pt_seed_manifest = Channel
        .of("dataset\tmodel\tdata_type\tcell_type\tchr\tcov_spec\tint_cov\tpt_seed\tk\tregion_file\tvariant_file\ttime_file")
        .concat(pt_seed_quasar.map { info, region_file, variant_file, time_file ->
            "${info.dataset}\t${info.model}\t${info.data_type}\t${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.int_cov}\t${info.pt_seed}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}"
        })
        .collectFile(
            name: "pseudotime-seed-quasar.tsv",
            newLine: true,
            sort: false,
            storeDir: "${projectDir}/output"
        )
    TABULATE_PT_SEED_EGENES(pt_seed_manifest, gene_map)
}
