include { BUILD_CYTOIMMGEN_SC_COVS } from '../modules/cytoimmgen'
include { RUN_QUASAR_SC as RUN_QUASAR_SC_CYTOIMMGEN } from '../modules/quasar'

// Main-effect p_glmm_sc per compartment x timepoint, and interactions on
// pooled (_all) compartments for timepoint, pseudotime and division state.
workflow CYTOIMMGEN_QTL {
    take:
    sc_counts
    pb_expr_covs
    geno_pcs
    harmonised_h5ad
    anno_bed
    bed_files

    main:
    def cyto_tps = params.dataset_configs.cytoimmgen.timepoints.collect { it.toString() }
    def cyto_timepoint_int = cyto_tps.drop(1).collect { "timepoint_${it}" }.join(",")
    def cyto_info_key = { info ->
        [info.dataset, info.cell_type, info.cell_frac, info.indiv_frac, info.min_n_cells, info.count_frac]
    }

    cyto_build_input = sc_counts
        .filter { _type, info, _pheno -> info.dataset == "cytoimmgen" && info.is_full_data }
        .map { _type, info, pheno -> tuple(info, pheno) }
        .combine(
            pb_expr_covs.filter { info, _covs -> info.dataset == "cytoimmgen" },
            by: 0
        )
        .map { info, pheno, pb_covs ->
            tuple(info.dataset, info, pheno, pb_covs)
        }
        .combine(geno_pcs, by: 0)
        .combine(harmonised_h5ad, by: 0)
        .flatMap { dataset, info, pheno, pb_covs, pcs, h5ad ->
            def specs
            if (info.cell_type.toString().endsWith("_all")) {
                specs = [
                    ["cyto+timepoint", "timepoint", cyto_timepoint_int],
                    ["cyto+pseudotime", "pseudotime", "pseudotime"],
                    ["cyto+division_state", "division_state", "division_state"],
                ]
            } else {
                specs = [["cyto", "none", "none"]]
            }
            specs.collect { cov_spec, int_cov, interaction_cov ->
                tuple(
                    info + [
                        cov_spec: cov_spec,
                        int_cov: int_cov,
                        interaction_cov: interaction_cov,
                        data_type: "counts",
                        model: "p_glmm_sc",
                        k: "none",
                    ],
                    pheno,
                    pb_covs,
                    pcs,
                    h5ad,
                    cov_spec,
                    cyto_tps
                )
            }
        }

    cyto_covs = BUILD_CYTOIMMGEN_SC_COVS(cyto_build_input)

    cyto_quasar_input = cyto_covs
        .map { info, cov_spec, covs ->
            tuple(cyto_info_key(info), info, cov_spec, covs)
        }
        .combine(
            sc_counts
                .filter { _type, info, _pheno -> info.dataset == "cytoimmgen" }
                .map { _type, info, pheno -> tuple(cyto_info_key(info), pheno) },
            by: 0
        )
        .map { _key, info, cov_spec, covs, pheno ->
            tuple(info.dataset, info, cov_spec, covs, pheno)
        }
        .combine(anno_bed, by: 0)
        .combine(bed_files.filter { dataset, _chr, _bed -> dataset == "cytoimmgen" }, by: 0)
        .map { dataset, info, cov_spec, covs, pheno, anno, chr, bed ->
            tuple(
                info + [chr: chr],
                pheno,
                covs,
                anno,
                bed,
                covs
            )
        }

    quasar = RUN_QUASAR_SC_CYTOIMMGEN(cyto_quasar_input)

    emit:
    quasar
}
