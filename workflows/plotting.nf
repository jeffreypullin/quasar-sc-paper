include {
    PLOT_POWER ; PLOT_POWER_QUASAR ; PLOT_POWER_COVS ; PLOT_CONVERGENCE ; PLOT_PERM ;
    PLOT_PERM_CELL_FRAC ; PLOT_PERM_COUNT_FRAC ; PLOT_PERM_PB ;
    PLOT_TIME ; PLOT_CONCORDANCE ; PLOT_UNIQUE_SC_EGENES ;
    PLOT_UNIQUE_PB_EGENES ; PLOT_EGENE_SIG_MODEL ; PLOT_POWER_CELL_FRAC ;
    PLOT_POWER_COUNT_FRAC ; PLOT_POWER_INDIV_FRAC ; PLOT_POWER_BOTH_FRAC ;
    PLOT_POWER_GENE_PROP ; PLOT_METHOD_SCATTER ; PLOT_GENE_PROPERTIES ;
    PLOT_PERM_INT ; PLOT_INT_OUTPUT ; PLOT_GWAS_OUTPUT ; PLOT_PERM_GLOBAL ;
    PLOT_PERM_GWAS ; PLOT_PVALUE_SCATTER ; PLOT_MAIN_VS_INT ;
    PLOT_INT_TIME ; PLOT_QUASAR_INT_TIME ; PLOT_SC_INT_OUTPUT ;
    PLOT_QUASAR_SC_INT_OUTPUT ; PLOT_PERM_SC_INT ; PLOT_PERM_SC_INT_WITHIN ;
    PLOT_PERM_SC_INT_LOGCOUNT_BINS ; PLOT_PERM_SC_GROUPED ;
    PLOT_GROUPED_SC_OUTPUT ; PLOT_GROUPED_VS_INT ; PLOT_GROUPED_INT_TIME ;
    PLOT_SC_PGLMM_VS_LMM ; PLOT_PERM_SC_INT_GLOBAL ; PLOT_SC_INT_FIGURES ;
    PLOT_UNIQUE_SC_EGENE_FIGURES ; PLOT_UNIQUE_PB_EGENE_FIGURES ;
    PLOT_METACELL_OUTPUT ; PLOT_CSAQTL_OUTPUT ; PLOT_PC_GWAS_OUTPUT ;
    PLOT_PC_GWAS_COMPARISON ; RUN_COLOC ; PLOT_CELLS_PER_INDIV ;
    PLOT_POWER_N_CELLS_FILTER ; PLOT_POWER_OFFSET ; PLOT_ZSCORE_SCATTER ;
    PLOT_TREMOR_TENSORQTL ; PLOT_DATASET_QC
} from '../modules/analysis'

workflow RUN_PLOTS {
    take:
    gene_properties_file
    pb_quasar_file
    sc_quasar_file
    sc_offset_quasar_file
    seacells_file
    sc_quasar_gwas_file
    pb_quasar_gwas_file
    csaqtl_quasar_file
    pc_gwas_file
    pc_sc_gwas_file
    sc_quasar_gwas_perm_file
    pb_quasar_gwas_perm_file
    perm_int_pb_quasar_file
    sc_quasar_perm_int_file
    sc_quasar_perm_int_within_file
    sc_quasar_perm_int_logcount_file
    perm_sc_quasar_file
    perm_pb_quasar_file
    saigeqtl_file
    cluster_sizes
    tremor_tensorqtl_compare
    dataset_qc
    active_datasets

    main:
    // PLOT_CLUMPED(clumped_pb_quasar_file, clumped_sc_quasar_file)
    // PLOT_POWER(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    PLOT_POWER_QUASAR(pb_quasar_file, sc_quasar_file)
    // PLOT_POWER_N_CELLS_FILTER(pb_quasar_file, sc_quasar_file)
    // PLOT_POWER_OFFSET(sc_offset_quasar_file)
    // unique_egenes_tsv = PLOT_UNIQUE_SC_EGENES(pb_quasar_file, sc_quasar_file)
    // PLOT_EGENE_SIG_MODEL(pb_quasar_file, sc_quasar_file)
    // PLOT_POWER_CELL_FRAC(pb_quasar_file, sc_quasar_file)
    // TEMP: disabled with joint downsampling
    // PLOT_POWER_BOTH_FRAC(pb_quasar_file, sc_quasar_file)
    // PLOT_POWER_COUNT_FRAC(pb_quasar_file, sc_quasar_file)
    // PLOT_POWER_INDIV_FRAC(pb_quasar_file, sc_quasar_file)
    // PLOT_POWER_GENE_PROP(pb_quasar_file, sc_quasar_file)
    // PLOT_METHOD_SCATTER(pb_quasar_file, sc_quasar_file)

    // unique_genotype_dosages = EXTRACT_GENOTYPES(
    //     all_bed.combine(unique_egenes_tsv)
    // )
    // unique_sc_logcounts = EXTRACT_SC_LOGCOUNTS(
    //     sc_preprocess_input
    //         .filter { it[0].is_full_data && it[0].cell_type in ["Plasma", "B_IN", "CD4_NC"] }
    //         .combine(unique_egenes_tsv)
    // )
    // unique_sc_egene_figures_input = unique_sc_logcounts
    //     .map { info, logcounts -> tuple(info.dataset, info, logcounts) }
    //     .combine(unique_genotype_dosages, by: 0)
    //     .combine(unique_egenes_tsv)
    //     .map { dataset, info, logcounts, geno, leads ->
    //         tuple(info, leads, logcounts, geno)
    //     }
    // PLOT_UNIQUE_SC_EGENE_FIGURES(unique_sc_egene_figures_input)

    // unique_pb_egenes_tsv = PLOT_UNIQUE_PB_EGENES(pb_quasar_file, sc_quasar_file)
    // unique_pb_genotype_dosages = EXTRACT_GENOTYPES_PB_UNIQUE(
    //     all_bed.combine(unique_pb_egenes_tsv)
    // )
    // unique_pb_sc_logcounts = EXTRACT_SC_LOGCOUNTS_PB_UNIQUE(
    //     sc_preprocess_input
    //         .filter { it[0].is_full_data && it[0].cell_type in ["Plasma", "B_IN", "CD4_NC"] }
    //         .combine(unique_pb_egenes_tsv)
    // )
    // unique_pb_egene_figures_input = unique_pb_sc_logcounts
    //     .map { info, logcounts -> tuple(info.dataset, info, logcounts) }
    //     .combine(unique_pb_genotype_dosages, by: 0)
    //     .combine(unique_pb_egenes_tsv)
    //     .map { dataset, info, logcounts, geno, leads ->
    //         tuple(info, leads, logcounts, geno)
    //     }
    // PLOT_UNIQUE_PB_EGENE_FIGURES(unique_pb_egene_figures_input)

    // PLOT_POWER_COVS(sc_quasar_file)
    // PLOT_CONVERGENCE(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    PLOT_PERM(perm_sc_quasar_file)
    // PLOT_PERM_CELL_FRAC(perm_sc_quasar_file, perm_pb_quasar_file)
    // PLOT_PERM_COUNT_FRAC(perm_sc_quasar_file, perm_pb_quasar_file)
    PLOT_PERM_PB(perm_pb_quasar_file)

    if ("tremor" in active_datasets) {
        tremor_plot_files = tremor_tensorqtl_compare
            .map { dataset, cell_type, sample, summary -> tuple(sample, summary) }
            .toList()
            .filter { it.size() > 0 }
            .map { pairs ->
                tuple(
                    pairs.collect { it[0] },
                    pairs.collect { it[1] }
                )
            }
        PLOT_TREMOR_TENSORQTL(
            tremor_plot_files.map { samples, summaries -> samples },
            tremor_plot_files.map { samples, summaries -> summaries }
        )
    }

    qc_files = dataset_qc
        .toList()
        .filter { it.size() > 0 }
        .map { rows ->
            tuple(
                rows.collect { it[1] },
                rows.collect { it[2] }
            )
        }
    PLOT_DATASET_QC(
        qc_files.map { counts, cells -> counts },
        qc_files.map { counts, cells -> cells }
    )

    // PLOT_PERM_GLOBAL(perm_pb_quasar_file)
    // PLOT_PERM_GWAS(sc_quasar_gwas_perm_file, pb_quasar_gwas_perm_file)
    // PLOT_PERM_INT(perm_int_pb_quasar_file)
    // PLOT_INT_OUTPUT(pb_quasar_file)
    // PLOT_SC_INT_OUTPUT(sc_quasar_file, castie_file)
    // PLOT_QUASAR_SC_INT_OUTPUT(sc_quasar_file)
    // PLOT_GROUPED_SC_OUTPUT(sc_quasar_file)
    // PLOT_GROUPED_INT_TIME(sc_quasar_file)
    // PLOT_GROUPED_VS_INT(sc_quasar_file, castie_file)
    // PLOT_SC_PGLMM_VS_LMM(sc_quasar_file)
    // PLOT_PERM_SC_GROUPED(perm_sc_quasar_file)
    // PLOT_PERM_SC_INT_GLOBAL(perm_sc_quasar_file)
    // PLOT_PERM_SC_INT(sc_quasar_perm_int_file)
    // PLOT_SC_INT_FIGURES(sc_int_figures_input)
    // PLOT_METACELL_OUTPUT(sc_quasar_file, seacells_file)
    // PLOT_CSAQTL_OUTPUT(csaqtl_quasar_file)
    // PLOT_PC_GWAS_OUTPUT(pc_gwas_file, pc_sc_gwas_file)
    // PLOT_PERM_SC_INT_WITHIN(sc_quasar_perm_int_within_file)
    // PLOT_PERM_SC_INT_LOGCOUNT_BINS(sc_quasar_perm_int_logcount_file)
    // PLOT_TIME(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    // PLOT_CELLS_PER_INDIV(
    //     cluster_sizes.map { dataset, sizes, cells_per_indiv -> cells_per_indiv }
    // )
    // PLOT_PC_GWAS_COMPARISON(pc_gwas_file, pc_sc_gwas_file)
    // PLOT_CONCORDANCE(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    // PLOT_GENE_PROPERTIES(gene_properties_file)
    // PLOT_GWAS_OUTPUT(sc_quasar_gwas_file, pb_quasar_gwas_file)
    // PLOT_PVALUE_SCATTER(pb_quasar_file, sc_quasar_file)
    // PLOT_ZSCORE_SCATTER(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    // PLOT_MAIN_VS_INT(sc_quasar_file)
    // PLOT_INT_TIME(sc_quasar_file, castie_file)
    // PLOT_QUASAR_INT_TIME(sc_quasar_file)
    // RUN_COLOC(sc_quasar_file)
}
