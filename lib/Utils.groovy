class Utils {

    static attachGeneProperties(ch, gene_properties) {
        ch
            .map { info, region, variant, time ->
                [[dataset: info.dataset, cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac, min_n_cells: info.min_n_cells, count_frac: info.count_frac],
                 info, region, variant, time]
            }
            .combine(
                gene_properties.map { info, gp ->
                    [[dataset: info.dataset, cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac, min_n_cells: info.min_n_cells, count_frac: info.count_frac], gp]
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

    // Donors left by apply_cell_selection (indiv_frac assumed 1): each donor keeps
    // max(1, floor(cell_frac * n)) cells, then must have >= min_n_cells.
    static int nDonorsAfterSelection(cells_per_indiv, info) {
        def label = info.cell_type.replace("_", " ")
        java.nio.file.Files.readAllLines(cells_per_indiv)
            .drop(1)
            .collect { it.split("\t") }
            .findAll { it[0] == label && (it[2] as int) > 0 }
            .count { row ->
                def n = row[2] as int
                def n_ds = info.cell_frac < 1.0d ? Math.max(1, Math.floor(info.cell_frac * n) as int) : n
                n_ds >= info.min_n_cells
            } as int
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
            && info.min_n_cells < 0
            && !isJointBothFrac(info)
    }

    // Proportional cell-fraction series on the three example cell types.
    static boolean isCellFracSeries(info) {
        info.dataset == "onek1k"
            && info.cell_type in ["Plasma", "B_IN", "CD4_NC"]
            && info.indiv_frac == 1.0d
            && info.count_frac == 1.0d
            && info.min_n_cells < 0
            && !isJointBothFrac(info)
    }

    // Cell-fraction × min_n_cells (keep donors with >= 5 or >= 10 cells after subsample).
    // Plasma skipped for now.
    static boolean isCellFracDonorFilterSeries(info) {
        info.dataset == "onek1k"
            && info.cell_type in ["B_IN", "CD4_NC"]
            && info.indiv_frac == 1.0d
            && info.count_frac == 1.0d
            && info.min_n_cells in [5, 10]
            && !isJointBothFrac(info)
    }

    // Count-depth series on the three main cell types, including the full-data point.
    static boolean isCountFracSeries(info) {
        info.dataset == "onek1k"
            && info.cell_type in ["Plasma", "B_IN", "CD4_NC"]
            && info.cell_frac == 1.0d
            && info.indiv_frac == 1.0d
            && info.min_n_cells < 0
            && !isJointBothFrac(info)
    }

    static String ablateSuffix(info) {
        info.ablate ? "-ablate-${info.ablate}".toString() : ""
    }

    // SPA is ablated by a cutoff no finite standardised score can exceed.
    static String ablateFlag(info) {
        if (!info.ablate) return ""
        info.ablate == "spa" ? " --spa-cutoff 1e300" : " --ablate-${info.ablate}".toString()
    }

    // Tremor: full-data lm. tmm_int uses the TensorQTL covariates;
    // logcounts stays on the OneK1K covariate set (sex, age, expr PC1–5, geno PC1–6).
    // OneK1K: full-data lm, lm_cw, nb_glm, and p_glmm_id on the main clusters.
    static boolean keepPbQuasarTemp(info) {
        if (info.int_cov != "none" || !info.is_full_data) {
            return false
        }
        if (info.dataset == "tremor") {
            return info.model == "lm" && info.pb_type in ["tmm_int", "logcounts"]
        }
        info.dataset == "onek1k"
            && !(info.cell_type in ["B_all", "T_all"])
            && ((info.model == "lm" && info.pb_type == "logcounts")
                || (info.model == "lm_cw" && info.pb_type == "logcounts")
                || (info.model == "nb_glm" && info.pb_type == "counts")
                || (info.model == "p_glmm_id" && info.pb_type == "counts"))
    }

    // Tremor full-data lm genotype permutations, plus the OneK1K
    // pseudobulk-lm cell-downsampling series used for the cell-fraction QQ plot.
    static boolean keepPbPermTemp(info) {
        if (info.dataset == "tremor") {
            return info.model == "lm"
                && info.pb_type in ["tmm_int", "logcounts"]
                && info.int_cov == "none"
                && info.is_full_data
        }
        info.model == "lm"
            && info.pb_type == "logcounts"
            && info.int_cov == "none"
            && isCellFracSeries(info)
    }

    // OneK1K uses the selection already applied to sc_quasar_input_full.
    // Tremor stays on full-data counts with p_glmm_sc and logcounts with lmm_sc,
    // each with bulk PCA and TensorQTL covariates.
    static boolean keepScQuasarTemp(info) {
        if (info.dataset == "onek1k") {
            return true
        }
        info.dataset == "tremor"
            && info.is_full_data
            && info.k == "none"
            && info.int_cov == "none"
            && info.cov_spec in ["bulk_pca", "tqtl"]
            && !info.ablate
            && ((info.model == "p_glmm_sc" && info.data_type == "counts")
                || (info.model == "lmm_sc" && info.data_type == "log_counts"))
    }

    // Genotype permutations for the non-zero-fraction QQ plot. Tremor stays on
    // full data. OneK1K also runs the p_glmm_sc cell-downsampling series used
    // for the cell-fraction QQ plot.
    static boolean keepScPermTemp(info) {
        def is_sc_perm_model =
            (info.model == "p_glmm_sc" && info.data_type == "counts")
                || (info.model == "lmm_sc" && info.data_type == "log_counts")
        if (info.dataset == "onek1k"
            && !(info.cell_type in ["B_all", "T_all"])
            && info.k == "none"
            && info.int_cov == "none"
            && info.cov_spec in ["bulk_pca", "tqtl"]
            && !info.ablate
            && is_sc_perm_model) {
            return info.is_full_data
                || (info.model == "p_glmm_sc"
                    && info.data_type == "counts"
                    && info.cov_spec == "bulk_pca"
                    && isCellFracSeries(info))
        }
        info.dataset == "tremor"
            && !(info.cell_type in ["B_all", "T_all"])
            && info.is_full_data
            && info.k == "none"
            && info.int_cov == "none"
            && info.cov_spec in ["bulk_pca", "tqtl"]
            && !info.ablate
            && is_sc_perm_model
    }

    // Genotype permutations of the OneK1K B_all pseudotime interaction,
    // including the K = 5 grouped fit.
    static boolean keepScIntPermTemp(info) {
        def is_sc_perm_model =
            (info.model == "p_glmm_sc" && info.data_type == "counts")
                || (info.model == "lmm_sc" && info.data_type == "log_counts")
        info.dataset == "onek1k"
            && info.cell_type == "B_all"
            && info.is_full_data
            && info.cov_spec == "bulk_pca+pseudotime"
            && info.int_cov == "pseudotime"
            && info.k in ["none", 5]
            && !info.ablate
            && is_sc_perm_model
    }

    // Which subsample x model x covariate combinations reach single-cell QuASAR.
    static boolean keepScQuasarCovSeries(info) {
        if (info.dataset == "tremor") return info.is_full_data
        info.is_full_data
           || (info.model in ["p_glmm_sc", "lmm_sc"]
               && info.cov_spec == "bulk_pca"
               && info.indiv_frac == 1.0d
               && info.count_frac == 1.0d)
           || (info.model in ["p_glmm_sc", "lmm_sc"]
               && info.cov_spec == "bulk_pca"
               && info.cell_type == "CD4_NC"
               && info.cell_frac == 1.0d
               && info.count_frac == 1.0d)
           || (info.model == "p_glmm_sc"
               && info.cov_spec == "bulk_pca"
               && isJointBothFrac(info))
           || (info.model in ["p_glmm_sc", "lmm_sc"]
               && info.cov_spec == "bulk_pca"
               && info.min_n_cells >= 0
               && info.count_frac == 1.0d)
           || (info.model in ["p_glmm_sc", "lmm_sc"]
               && info.cov_spec == "bulk_pca"
               && info.cell_frac == 1.0d
               && info.indiv_frac == 1.0d
               && info.min_n_cells < 0
               && info.count_frac != 1.0d)
    }

    // Tremor/Randolph SC mapping use configured cell types.
    // OneK1K main-effect SC is all clusters except B_all / T_all,
    // plus the count-frac and cell-fraction series on Plasma / B_IN / CD4_NC,
    // plus CD4_NC donor-fraction subsampling and the B_IN / CD4_NC
    // cell-fraction grid with a 5- or 10-cell donor filter.
    static boolean keepScQuasarCellType(info) {
        (info.dataset == "tremor" && info.is_full_data)
           || isRandolphDataset(info.dataset)
           || info.cell_type in ["B_all", "T_all"]
           || (info.dataset == "onek1k"
               && info.model in ["p_glmm_sc", "lmm_sc"]
               && !(info.cell_type in ["B_all", "T_all"])
               && info.is_full_data
               && info.k == "none"
               && info.int_cov == "none")
           || (info.model in ["p_glmm_sc", "lmm_sc"]
               && info.cell_type in ["Plasma", "B_IN", "CD4_NC"]
               && info.k == "none"
               && info.int_cov == "none"
               && info.cell_frac == 1.0d
               && info.indiv_frac == 1.0d
               && info.min_n_cells < 0)
           || (info.model in ["p_glmm_sc", "lmm_sc"]
               && info.cell_type in ["Plasma", "B_IN", "CD4_NC"]
               && info.k == "none"
               && info.int_cov == "none"
               && info.cell_frac == 1.0d
               && info.indiv_frac == 1.0d
               && info.min_n_cells >= 0
               && info.count_frac == 1.0d)
           || (info.model in ["p_glmm_sc", "lmm_sc"]
               && info.cell_type == "CD4_NC"
               && info.k == "none"
               && info.int_cov == "none"
               && info.count_frac == 1.0d
               && info.min_n_cells < 0
               && !isJointBothFrac(info)
               && (info.cell_frac != 1.0d || info.indiv_frac != 1.0d))
           || (info.dataset == "onek1k"
               && info.model == "p_glmm_sc"
               && !(info.cell_type in ["B_all", "T_all"])
               && info.k == "none"
               && info.int_cov == "none"
               && info.indiv_frac == 1.0d
               && info.count_frac == 1.0d
               && info.min_n_cells < 0
               && !isJointBothFrac(info))
           || (info.model == "p_glmm_sc"
               && info.k == "none"
               && info.int_cov == "none"
               && isCellFracDonorFilterSeries(info))
    }

    // Full data, plus the B_all grouped interaction and the subsample series.
    static boolean keepScQuasarScale(info) {
        if (info.dataset == "tremor") return info.is_full_data
        info.is_full_data
           || (info.model == "p_glmm_sc"
               && info.k == 5
               && info.int_cov != "none"
               && info.cell_type in ["B_all"])
           || (info.model in ["p_glmm_sc", "lmm_sc"]
               && info.k == "none"
               && info.int_cov == "none"
               && info.indiv_frac == 1.0d
               && info.count_frac == 1.0d)
           || (info.model in ["p_glmm_sc", "lmm_sc"]
               && info.k == "none"
               && info.int_cov == "none"
               && info.cell_type == "CD4_NC"
               && info.cell_frac == 1.0d
               && info.count_frac == 1.0d)
           || (info.model == "p_glmm_sc"
               && info.k == "none"
               && info.int_cov == "none"
               && isJointBothFrac(info))
           || (info.model in ["p_glmm_sc", "lmm_sc"]
               && info.k == "none"
               && info.int_cov == "none"
               && info.min_n_cells >= 0
               && info.count_frac == 1.0d)
           || (info.model in ["p_glmm_sc", "lmm_sc"]
               && info.k == "none"
               && info.int_cov == "none"
               && info.cell_frac == 1.0d
               && info.indiv_frac == 1.0d
               && info.min_n_cells < 0
               && info.count_frac != 1.0d)
    }
}
