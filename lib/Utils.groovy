class Utils {

    static attachGeneProperties(ch, gene_properties) {
        ch
            .map { info, region, variant, time ->
                [[dataset: info.dataset, cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac, n_cells_target: info.n_cells_target, count_frac: info.count_frac],
                 info, region, variant, time]
            }
            .combine(
                gene_properties.map { info, gp ->
                    [[dataset: info.dataset, cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac, n_cells_target: info.n_cells_target, count_frac: info.count_frac], gp]
                },
                by: 0
            )
            .map { base, info, region, variant, time, gp -> [info, region, variant, time, gp] }
    }

    static combineWithPrunedSnps(ch, pruned_snps) {
        ch
            .map { info, region_file, variant_file, time_file ->
                [info.dataset, info.chr, info, region_file, variant_file, time_file]
            }
            .combine(pruned_snps, by: [0, 1])
            .map { _dataset, _chr, info, region_file, variant_file, time_file, prune_in ->
                [info, region_file, variant_file, time_file, prune_in]
            }
    }

    static boolean isRandolphDataset(ds) {
        ds.toString().startsWith("randolph_")
    }

    static boolean isJointBothFrac(info) {
        info.joint_both_frac == true
    }

    // CD4_NC donor-fraction series, including the full-data point.
    static boolean isDonorSubsampleSeries(info) {
        info.dataset == "onek1k"
            && info.cell_type == "CD4_NC"
            && info.cell_frac == 1.0d
            && info.count_frac == 1.0d
            && info.n_cells_target < 0
            && !isJointBothFrac(info)
    }

    // Proportional cell-fraction series on the three main cell types, including the full-data point.
    static boolean isCellFracSeries(info) {
        info.dataset == "onek1k"
            && info.cell_type in ["Plasma", "B_IN", "CD4_NC"]
            && info.indiv_frac == 1.0d
            && info.count_frac == 1.0d
            && info.n_cells_target < 0
            && !isJointBothFrac(info)
    }

    // Count-depth series on the three main cell types, including the full-data point.
    static boolean isCountFracSeries(info) {
        info.dataset == "onek1k"
            && info.cell_type in ["Plasma", "B_IN", "CD4_NC"]
            && info.cell_frac == 1.0d
            && info.indiv_frac == 1.0d
            && info.n_cells_target < 0
            && !isJointBothFrac(info)
    }

    static String ablateSuffix(info) {
        info.ablate ? "-ablate-${info.ablate}".toString() : ""
    }

    // TEMP: OneK1K / Randolph full-data PB lm / lm_cw / nb_glm / p_glmm_id for the power plot,
    // plus CD4_NC donor-fraction lm.
    static boolean keepPbQuasarTemp(info) {
        ((info.dataset == "onek1k" || isRandolphDataset(info.dataset))
            && info.model in ["lm", "lm_cw", "nb_glm", "p_glmm_id"]
            && !(info.cell_type in ["B_all", "T_all"])
            && info.is_full_data
            && info.int_cov == "none")
            || (info.model == "lm"
                && info.int_cov == "none"
                && isDonorSubsampleSeries(info))
    }

    // TEMP: OneK1K full-data nb_glm genotype perms on all cell types,
    // p_glmm_id on the three main cell types, plus lm on the CD4_NC donor-fraction
    // and three-cell-type cell-fraction series.
    static boolean keepPbPermTemp(info) {
        (info.dataset == "onek1k"
            && (info.model == "nb_glm"
                || (info.model == "p_glmm_id" && info.cell_type in ["Plasma", "B_IN", "CD4_NC"]))
            && !(info.cell_type in ["B_all", "T_all"])
            && info.is_full_data
            && info.int_cov == "none")
            || (info.model == "lm"
                && info.int_cov == "none"
                && (isDonorSubsampleSeries(info) || isCellFracSeries(info) || isCountFracSeries(info)))
    }

    // TEMP: OneK1K / Randolph full-data single-cell P-GLMM for the QuASAR power plot,
    // plus the CD4_NC donor-fraction series, plus the cell-level covariate
    // comparison on Plasma / B_IN / CD4_NC, plus the full-data B_all pseudotime interaction
    // (continuous and K = 5 binned) and its matching main-effect run for timing comparison.
    static boolean keepScQuasarTemp(info) {
        ((info.dataset == "onek1k" || isRandolphDataset(info.dataset))
            && info.model == "p_glmm_sc"
            && info.data_type == "counts"
            && info.k == "none"
            && info.int_cov == "none"
            && !info.ablate
            && ((info.cov_spec == "bulk_pca"
                    && ((!(info.cell_type in ["B_all", "T_all"]) && info.is_full_data)
                        || isDonorSubsampleSeries(info)))
                || (info.cov_spec in ["sc_pca", "mean_sc_pca", "bulk_pca+pct_mito", "bulk_pca+cell_cycle"]
                    && info.cell_type in ["Plasma", "B_IN", "CD4_NC"]
                    && info.is_full_data)))
            || (info.dataset == "onek1k"
                && info.model == "p_glmm_sc"
                && info.data_type in ["counts", "sct_counts"]
                && info.k in ["none", 5]
                && info.int_cov == "pseudotime"
                && info.cell_type == "B_all"
                && info.is_full_data
                && !info.ablate)
            || (info.dataset == "onek1k"
                && info.model == "p_glmm_sc"
                && info.data_type == "counts"
                && info.k == "none"
                && info.int_cov == "none"
                && info.cov_spec == "bulk_pca"
                && info.cell_type == "B_all"
                && info.is_full_data
                && !info.ablate)
    }

    // TEMP: OneK1K / Randolph p_glmm_sc genotype perms on full data for all cell types,
    // plus the OneK1K cell-fraction series and CD4_NC donor-fraction, for calibration plots.
    static boolean keepScPermTemp(info) {
        ((info.dataset == "onek1k" || isRandolphDataset(info.dataset))
            && info.model == "p_glmm_sc"
            && info.data_type == "counts"
            && info.cov_spec == "bulk_pca"
            && info.k == "none"
            && info.int_cov == "none"
            && !info.ablate
            && ((!(info.cell_type in ["B_all", "T_all"]) && info.is_full_data)
                || isCellFracSeries(info)
                || isDonorSubsampleSeries(info)
                || isCountFracSeries(info)))
    }
}
