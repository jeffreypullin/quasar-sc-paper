include {
    PLOT_POWER ; PLOT_POWER_QUASAR ; PLOT_POWER_COVS ; PLOT_CONVERGENCE ; PLOT_PERM_SC ;
    PLOT_PERM_CELL_FRAC ; PLOT_PERM_CELL_FRAC_MIN_CELLS ; PLOT_PERM_CELL_FRAC_DONOR_SPARSITY ;
    PLOT_PERM_COUNT_FRAC ;
    PLOT_PERM_DONOR_SUBSAMPLE ; PLOT_PERM_PB ;
    PLOT_TIME ; PLOT_TIME_QUASAR ; PLOT_CONCORDANCE ; PLOT_UNIQUE_SC_EGENES ;
    PLOT_UNIQUE_PB_EGENES ; PLOT_EGENE_SIG_MODEL ; PLOT_POWER_CELL_FRAC ;
    PLOT_CELL_SUBSAMPLE_TIME ; PLOT_DONOR_SUBSAMPLE_TIME ; PLOT_TIME_NCELLS_SCATTER ;
    PLOT_POWER_COUNT_FRAC ; PLOT_POWER_INDIV_FRAC ; PLOT_POWER_BOTH_FRAC ;
    PLOT_POWER_GENE_PROP ; PLOT_METHOD_SCATTER ; PLOT_GENE_PROPERTIES ;
    PLOT_PERM_INT ; PLOT_INT_OUTPUT ; PLOT_GWAS_OUTPUT ; PLOT_PERM_GLOBAL ;
    PLOT_PERM_GWAS ; PLOT_PVALUE_SCATTER ; PLOT_MAIN_VS_INT ;
    PLOT_INT_TIME ; PLOT_QUASAR_INT_TIME ; PLOT_SC_INT_OUTPUT ;
    PLOT_QUASAR_SC_INT_OUTPUT ; PLOT_PERM_SC_INT ; PLOT_PERM_SC_INT_WITHIN ;
    PLOT_PERM_SC_INT_LOGCOUNT_BINS ; PLOT_PERM_SC_ALL ; PLOT_PERM_SC_GROUPED ;
    PLOT_GROUPED_SC_OUTPUT ; PLOT_GROUPED_VS_INT ; PLOT_GROUPED_INT_TIME ;
    PLOT_SC_PGLMM_VS_LMM ; PLOT_PERM_SC_INT_GLOBAL ; PLOT_SC_INT_FIGURES ;
    PLOT_UNIQUE_SC_EGENE_FIGURES ; PLOT_UNIQUE_PB_EGENE_FIGURES ;
    PLOT_METACELL_OUTPUT ; PLOT_CSAQTL_OUTPUT ; PLOT_PC_GWAS_OUTPUT ;
    PLOT_PC_GWAS_COMPARISON ; RUN_COLOC ; PLOT_CELLS_PER_INDIV ;
    PLOT_POWER_OFFSET ; PLOT_ZSCORE_SCATTER ;
    PLOT_TREMOR_TENSORQTL ; PLOT_TREMOR_POWER_QUASAR ; PLOT_TREMOR_PERM_NZ_FRAC ;
    PLOT_DATASET_QC ; PLOT_MAF ; PLOT_ABLATE_DONOR_TOTALS ; PLOT_PROP_N_CELLS ;
    PLOT_GENE_FILTER_INTERSECTION
} from '../modules/analysis'
include { EXTRACT_GENOTYPES ; EXTRACT_SC_LOGCOUNTS } from '../modules/preprocess'

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
    skipped_combos_file
    saigeqtl_file
    castie_file
    cluster_sizes
    tremor_tensorqtl_compare
    dataset_qc
    active_datasets
    gene_map
    all_bed
    sc_preprocess_input
    sc_covs
    gene_filter_stats

    main:
    // PLOT_CLUMPED(clumped_pb_quasar_file, clumped_sc_quasar_file)
    // PLOT_POWER(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    //PLOT_ABLATE_DONOR_TOTALS(sc_quasar_file)
    // PLOT_POWER_QUASAR(pb_quasar_file, sc_quasar_file)
    // PLOT_CELLS_PER_INDIV(
    //     cluster_sizes.map { _dataset, _sizes, cells_per_indiv -> cells_per_indiv }
    // )
    // PLOT_TIME_QUASAR(pb_quasar_file, sc_quasar_file)
    // PLOT_MAF(pb_quasar_file, sc_quasar_file)
    // PLOT_CELL_SUBSAMPLE_TIME(pb_quasar_file, sc_quasar_file)
    // PLOT_DONOR_SUBSAMPLE_TIME(pb_quasar_file, sc_quasar_file)
    // PLOT_TIME_NCELLS_SCATTER(
    //     pb_quasar_file,
    //     sc_quasar_file,
    //     cluster_sizes.map { _dataset, sizes, _cells_per_indiv -> sizes }
    // )
    //PLOT_POWER_OFFSET(sc_offset_quasar_file)
    // unique_egenes_tsv = PLOT_UNIQUE_SC_EGENES(pb_quasar_file, sc_quasar_file)
    // PLOT_EGENE_SIG_MODEL(pb_quasar_file, sc_quasar_file)
    PLOT_POWER_CELL_FRAC(pb_quasar_file, sc_quasar_file)
    // TEMP: disabled with joint downsampling
    // PLOT_POWER_BOTH_FRAC(pb_quasar_file, sc_quasar_file)
    //PLOT_POWER_COUNT_FRAC(pb_quasar_file, sc_quasar_file)
    PLOT_POWER_INDIV_FRAC(pb_quasar_file, sc_quasar_file)
    gene_filter_stats
        .multiMap { _dataset, stats, donor ->
            stats: stats
            donor: donor
        }
        .set { gene_filter_parts }
    PLOT_GENE_FILTER_INTERSECTION(gene_filter_parts.stats, gene_filter_parts.donor)
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
    // PLOT_PERM_SC(perm_sc_quasar_file)
    PLOT_PERM_CELL_FRAC(perm_sc_quasar_file, perm_pb_quasar_file)
    PLOT_PERM_CELL_FRAC_MIN_CELLS(perm_sc_quasar_file, skipped_combos_file)
    PLOT_PERM_CELL_FRAC_DONOR_SPARSITY(perm_sc_quasar_file)
    // PLOT_PERM_COUNT_FRAC(perm_sc_quasar_file, perm_pb_quasar_file)
    cluster_sizes
        .map { _dataset, _sizes, cells_per_indiv -> cells_per_indiv }
        .multiMap { cells_per_indiv ->
            donor_subsample: cells_per_indiv
            prop_n_cells: cells_per_indiv
        }
        .set { cells_per_indiv_ch }
    PLOT_PERM_DONOR_SUBSAMPLE(
         perm_sc_quasar_file,
         perm_pb_quasar_file,
         sc_quasar_file,
         cells_per_indiv_ch.donor_subsample
    )
    // PLOT_PERM_PB(perm_pb_quasar_file)
    PLOT_PROP_N_CELLS(
        perm_sc_quasar_file,
        cells_per_indiv_ch.prop_n_cells
    )

    // if ("tremor" in active_datasets) {
    //     tremor_plot_files = tremor_tensorqtl_compare
    //         .map { dataset, cell_type, sample, summary -> tuple(sample, summary) }
    //         .toList()
    //         .filter { it.size() > 0 }
    //         .map { pairs ->
    //             tuple(
    //                 pairs.collect { it[0] },
    //                 pairs.collect { it[1] }
    //             )
    //         }
    //     PLOT_TREMOR_TENSORQTL(
    //         tremor_plot_files.map { samples, summaries -> samples },
    //         tremor_plot_files.map { samples, summaries -> summaries }
    //     )
    //     PLOT_TREMOR_POWER_QUASAR(pb_quasar_file, sc_quasar_file)
    //     PLOT_TREMOR_PERM_NZ_FRAC(perm_sc_quasar_file, perm_pb_quasar_file)
    // }

    // qc_files = dataset_qc
    //     .toList()
    //     .filter { it.size() > 0 }
    //     .map { rows ->
    //         tuple(
    //             rows.collect { it[1] },
    //             rows.collect { it[2] },
    //             rows.collect { it[3] }
    //         )
    //     }
    // PLOT_DATASET_QC(
    //     qc_files.map { counts, cells, zeros -> counts },
    //     qc_files.map { counts, cells, zeros -> cells },
    //     qc_files.map { counts, cells, zeros -> zeros }
    // )

    // PLOT_PERM_GLOBAL(perm_pb_quasar_file)
    // PLOT_PERM_GWAS(sc_quasar_gwas_perm_file, pb_quasar_gwas_perm_file)
    // PLOT_PERM_INT(perm_int_pb_quasar_file)
    // PLOT_INT_OUTPUT(pb_quasar_file)
    // PLOT_SC_INT_OUTPUT(sc_quasar_file, castie_file, gene_map)
    // PLOT_QUASAR_SC_INT_OUTPUT(sc_quasar_file, gene_map)
    // PLOT_GROUPED_SC_OUTPUT(sc_quasar_file)
    // PLOT_GROUPED_INT_TIME(sc_quasar_file)
    // PLOT_GROUPED_VS_INT(sc_quasar_file, castie_file)
    // PLOT_SC_PGLMM_VS_LMM(sc_quasar_file)
    // PLOT_PERM_SC_GROUPED(perm_sc_quasar_file)
    // PLOT_PERM_SC_INT_GLOBAL(perm_sc_quasar_file)
    // PLOT_PERM_SC_INT(sc_quasar_perm_int_file)
    // Interaction eGene examples to plot: one row per gene / lead SNP.
    // sc_int_examples_tsv = file("$projectDir/resources/sc-int-examples.tsv")
    // sc_int_examples = channel.fromPath(sc_int_examples_tsv).splitCsv(header: true, sep: "\t")
    // sc_int_genotypes = EXTRACT_GENOTYPES(
    //     sc_int_examples
    //         .map { ex -> tuple(ex.dataset) }
    //         .unique()
    //         .combine(all_bed, by: 0)
    //         .map { dataset, bed -> tuple(dataset, bed, sc_int_examples_tsv) }
    // )
    // sc_int_logcounts = EXTRACT_SC_LOGCOUNTS(
    //     sc_preprocess_input
    //         .filter { info, _scf, _gene_props -> info.is_full_data }
    //         .map { info, scf, gene_props -> tuple([info.dataset, info.cell_type], info, scf, gene_props) }
    //         .combine(
    //             sc_int_examples.map { ex -> [[ex.dataset, ex.cell_type]] }.unique(),
    //             by: 0
    //         )
    //         .map { _key, info, scf, gene_props -> tuple(info, scf, gene_props, sc_int_examples_tsv) }
    // )
    // sc_int_figures_input = sc_int_examples
    //     .map { ex -> tuple([ex.dataset, ex.cell_type, "bulk_pca+${ex.int_cov}".toString()], ex) }
    //     .combine(
    //         sc_covs
    //             .filter { info, _cov_spec, _covs -> info.is_full_data }
    //             .map { info, cov_spec, covs -> tuple([info.dataset, info.cell_type, cov_spec], info, covs) },
    //         by: 0
    //     )
    //     .map { _key, ex, info, covs -> tuple([ex.dataset, ex.cell_type], info, ex, covs) }
    //     .combine(
    //         sc_int_logcounts.map { info, logcounts -> tuple([info.dataset, info.cell_type], logcounts) },
    //         by: 0
    //     )
    //     .map { _key, info, ex, covs, logcounts -> tuple(ex.dataset, info, ex, covs, logcounts) }
    //     .combine(sc_int_genotypes, by: 0)
    //     .map { _dataset, info, ex, covs, logcounts, geno -> tuple(info, ex, covs, logcounts, geno) }
    // PLOT_SC_INT_FIGURES(sc_int_figures_input)
    // PLOT_METACELL_OUTPUT(sc_quasar_file, seacells_file)
    // PLOT_CSAQTL_OUTPUT(csaqtl_quasar_file)
    // PLOT_PC_GWAS_OUTPUT(pc_gwas_file, pc_sc_gwas_file)
    // PLOT_PERM_SC_INT_WITHIN(sc_quasar_perm_int_within_file)
    // PLOT_PERM_SC_INT_LOGCOUNT_BINS(sc_quasar_perm_int_logcount_file)
    // PLOT_PERM_SC_ALL(
    //     PLOT_PERM_SC_INT_GLOBAL.out.qq
    //         .mix(PLOT_PERM_SC_INT.out.qq, PLOT_PERM_SC_INT_WITHIN.out.qq)
    //         .collect()
    // )
    // PLOT_TIME(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    // PLOT_PC_GWAS_COMPARISON(pc_gwas_file, pc_sc_gwas_file)
    // PLOT_CONCORDANCE(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    PLOT_GENE_PROPERTIES(gene_properties_file)
    // PLOT_GWAS_OUTPUT(sc_quasar_gwas_file, pb_quasar_gwas_file)
    // PLOT_PVALUE_SCATTER(pb_quasar_file, sc_quasar_file)
    // PLOT_ZSCORE_SCATTER(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    // PLOT_MAIN_VS_INT(sc_quasar_file)
    // PLOT_INT_TIME(sc_quasar_file, castie_file)
    // PLOT_QUASAR_INT_TIME(sc_quasar_file)
    // RUN_COLOC(sc_quasar_file)
}
