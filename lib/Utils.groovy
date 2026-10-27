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

    // Proportional cell-fraction series (Plasma skipped for now).
    static boolean isCellFracSeries(info) {
        info.dataset == "onek1k"
            && info.cell_type in ["B_IN", "CD4_NC"]
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

    // TEMP: donor-fraction + three-type cell-fraction lm.
    static boolean keepPbQuasarTemp(info) {
        info.model == "lm"
            && info.int_cov == "none"
            && (isDonorSubsampleSeries(info) || isCellFracSeries(info))
    }

    // TEMP: three-type cell-fraction lm genotype perms for cell-frac QQ plots.
    static boolean keepPbPermTemp(info) {
        info.model == "lm"
            && info.int_cov == "none"
            && isCellFracSeries(info)
    }

    // TEMP: donor-fraction + cell-fraction (+ donor-filtered cell-fraction) p_glmm_sc.
    static boolean keepScQuasarTemp(info) {
        info.model == "p_glmm_sc"
            && info.data_type == "counts"
            && info.cov_spec == "bulk_pca"
            && info.k == "none"
            && info.int_cov == "none"
            && !info.ablate
            && (isDonorSubsampleSeries(info)
                || isCellFracSeries(info)
                || isCellFracDonorFilterSeries(info))
    }

    // TEMP: same series for p_glmm_sc genotype perms.
    static boolean keepScPermTemp(info) {
        info.model == "p_glmm_sc"
            && info.data_type == "counts"
            && info.cov_spec == "bulk_pca"
            && info.k == "none"
            && info.int_cov == "none"
            && !info.ablate
            && (isDonorSubsampleSeries(info)
                || isCellFracSeries(info)
                || isCellFracDonorFilterSeries(info))
    }
}
