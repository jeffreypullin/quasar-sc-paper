#!/usr/bin/env nextflow

include { 
    EXTRACT_INDIV_IDS ; FILTER_VCF ; CONVERT_VCF_TO_BED ; 
    PRUNE_SNPS ; CONCAT_BED_FILES ; EXTRACT_GENOTYPES ;
    EXTRACT_GENOTYPES as EXTRACT_GENOTYPES_PB_UNIQUE ;
    COMPUTE_SC_COUNTS ; COMPUTE_SC_SCT_COUNTS ; COMPUTE_SC_LOGCOUNTS ; EXTRACT_SC_LOGCOUNTS ;
    EXTRACT_SC_LOGCOUNTS as EXTRACT_SC_LOGCOUNTS_PB_UNIQUE ;
    COMPUTE_PB_COUNTS ; 
    COMPUTE_CLUSTER_SIZES ; COMPUTE_PB_EXPR_COVS ; COMPUTE_SC_EXPR_COVS ; 
    COMPUTE_GENOTYPE_PCS ; CREATE_ANNOT_BED ; COMPUTE_PB_LOGCOUNTS ; 
    CREATE_GRM ; PREPARE_SLINGSHOT_ADATA ; PREPARE_SEACELLS_ADATA ; 
    RUN_SLINGSHOT ; JOIN_SC_INT_COVS ; JOIN_PB_INT_COVS ;
    DOWNLOAD_STARCAT_REF ; COMPUTE_STARCAT_COVS ; COMPUTE_CSAQTL_COUNTS ;
    COMPUTE_SC_PC_PHENO ; COMPUTE_SC_PC_SC_PHENO ; COMPUTE_SC_OFFSETS
} from './modules/preprocess'
include { COLLATE_SAIGEQTL_INPUT ; SAIGE_SUBSET_BED ; EXTRACT_GENES ; 
          EXTRACT_GENES as EXTRACT_CASTIE_GENES ;
          RUN_SAIGEQTL ; CONCAT_SAIGEQTL_TSVS ; COMPUTE_SAIGEQTL_POWER } from './modules/saige-qtl'
include { COLLATE_CASTIE_INPUT ; RUN_CASTIE ; CONCAT_CASTIE_TSVS } from './modules/castie'
include { COMBINE_PB_COVS ; COMBINE_SC_COVS ; FILTER_COVS ; 
          ANNOTATE_PHENO ; RUN_QUASAR_PB ; RUN_QUASAR_SC ; 
          RUN_QUASAR_SC_GWAS ; RUN_QUASAR_PB_GWAS ; RUN_QUASAR_SC_OFFSET ; 
          RUN_QUASAR_GWAS as RUN_CSAQTL_GWAS ;
          RUN_QUASAR_GWAS as RUN_PC_GWAS ;
          RUN_QUASAR_SC_PC_GWAS as RUN_PC_SC_GWAS } from './modules/quasar'
include { COMPUTE_CELL_GROUPS ; COMPUTE_SEACELLS ;
          COMPUTE_CELL_GROUPS as COMPUTE_CSAQTL_GROUPS } from './modules/cellgroups'
include { FILTER_VARIANTS ;
          FILTER_VARIANTS as FILTER_VARIANTS_SC ;
          FILTER_VARIANTS as FILTER_VARIANTS_PB ;
          FILTER_VARIANTS as FILTER_VARIANTS_PB_INT ;
          FILTER_VARIANTS as FILTER_VARIANTS_SC_INT ;
          FILTER_VARIANTS as FILTER_VARIANTS_SC_INT_WITHIN ;
          FILTER_VARIANTS as FILTER_VARIANTS_SC_INT_LOGCOUNT } from './modules/postprocess'
include { COMPUTE_QUASAR_POWER as COMPUTE_QUASAR_POWER_SC } from './modules/quasar'
include { COMPUTE_QUASAR_POWER as COMPUTE_QUASAR_POWER_PB } from './modules/quasar'
include { COMPUTE_QUASAR_POWER as COMPUTE_QUASAR_POWER_SC_OFFSET } from './modules/quasar'
include { RUN_QUASAR_PB as RUN_QUASAR_PB_PERM } from './modules/quasar'
include { RUN_QUASAR_SC as RUN_QUASAR_SC_PERM } from './modules/quasar'
include { RUN_QUASAR_PB as RUN_QUASAR_PB_PERM_INT } from './modules/quasar'
include { RUN_QUASAR_SC as RUN_QUASAR_SC_PERM_INT } from './modules/quasar'
include { RUN_QUASAR_SC as RUN_QUASAR_SC_PERM_WITHIN_INT } from './modules/quasar'
include { RUN_QUASAR_SC as RUN_QUASAR_SC_PERM_LOGCOUNT_INT } from './modules/quasar'
include { RUN_QUASAR_PB_GWAS as RUN_QUASAR_PB_GWAS_PERM } from './modules/quasar'
include { RUN_QUASAR_SC_GWAS as RUN_QUASAR_SC_GWAS_PERM } from './modules/quasar'
include { PERMUTE_BED ; PERMUTE_COV ; PERMUTE_SC_COV_WITHIN ;
          PERMUTE_SC_COV_LOGCOUNT_BINS ;
          PERMUTE_BED as PERMUTE_BED_GWAS ;
          PERMUTE_BED as PERMUTE_BED_PB ;
          PERMUTE_COV as PERMUTE_SC_COV } from './modules/permute'
include { 
    PLOT_POWER ; PLOT_POWER_COVS ; PLOT_CONVERGENCE ; PLOT_PERM ; 
    PLOT_PERM_CELL_FRAC ; PLOT_PERM_COUNT_FRAC ;
    PLOT_PERM_PB ; PLOT_TIME ; PLOT_CONCORDANCE ; PLOT_UNIQUE_SC_EGENES ;
    PLOT_UNIQUE_PB_EGENES ; PLOT_EGENE_SIG_MODEL ;     PLOT_POWER_CELL_FRAC ; PLOT_POWER_COUNT_FRAC ;
    PLOT_POWER_INDIV_FRAC ; PLOT_POWER_BOTH_FRAC ; PLOT_POWER_GENE_PROP ; PLOT_METHOD_SCATTER ;
    COMPUTE_GENE_PROPERTIES ; PLOT_GENE_PROPERTIES ; PLOT_PERM_INT ;
    PLOT_INT_OUTPUT ; PLOT_GWAS_OUTPUT ; CREATE_EXAMPLE_DATA ;
    PLOT_PERM_GLOBAL ; PLOT_PERM_GWAS ; PLOT_PVALUE_SCATTER ; PLOT_MAIN_VS_INT ;
    PLOT_INT_TIME ; PLOT_QUASAR_INT_TIME ;
    PLOT_SC_INT_OUTPUT ; PLOT_QUASAR_SC_INT_OUTPUT ; PLOT_PERM_SC_INT ; PLOT_PERM_SC_INT_WITHIN ;
    PLOT_PERM_SC_INT_LOGCOUNT_BINS ; 
    // CLUMP_VARIANTS ; 
    PLOT_PERM_SC_GROUPED ; PLOT_GROUPED_SC_OUTPUT ; PLOT_GROUPED_VS_INT ;
    PLOT_GROUPED_INT_TIME ;
    PLOT_SC_PGLMM_VS_LMM ; 
    PLOT_PERM_SC_INT_GLOBAL ; // PLOT_CLUMPED ; 
    PLOT_SC_INT_FIGURES ; PLOT_UNIQUE_SC_EGENE_FIGURES ;
    PLOT_UNIQUE_PB_EGENE_FIGURES ;
    PLOT_METACELL_OUTPUT ; PLOT_CSAQTL_OUTPUT ;     PLOT_PC_GWAS_OUTPUT ;
    PLOT_PC_GWAS_COMPARISON ; RUN_COLOC ; PLOT_CELLS_PER_INDIV ; PLOT_POWER_N_CELLS_FILTER ;
    PLOT_POWER_OFFSET ; PLOT_ZSCORE_SCATTER
    } from './modules/analysis'
// include { CLUMP_VARIANTS as CLUMP_VARIANTS_PB } from "./modules/analysis"
// include { CLUMP_VARIANTS as CLUMP_VARIANTS_SC } from "./modules/analysis"
include { COMPUTE_CONVERGENCE as COMPUTE_CONVERGENCE_PB } from './modules/analysis'
include { COMPUTE_CONVERGENCE as COMPUTE_CONVERGENCE_SC } from './modules/analysis'

workflow {

    def active_datasets = (params.datasets instanceof List)
        ? params.datasets.collect { it.toString() }
        : params.datasets.toString().split(',').collect { it.trim() }.findAll { it }

    active_datasets.each { name ->
        if (!params.dataset_configs.containsKey(name)) {
            error "Unknown dataset '${name}'. Available: ${params.dataset_configs.keySet().join(', ')}"
        }
    }
    def active_configs = params.dataset_configs.findAll { name, dcfg -> name in active_datasets }

    // Joint cell × individual fraction grid for CD4_NC power heatmap.
    def joint_fracs = [0.1d, 0.25d, 0.5d, 0.75d, 1.0d]
    def is_joint_both_frac = { info -> info.joint_both_frac == true }
    def is_starcat_int_cov = { int_cov -> int_cov != null && int_cov.toString().startsWith("starcat_") }
    def is_int_sc_cov = { int_cov ->
        int_cov == "none" || int_cov == "pseudotime" || is_starcat_int_cov(int_cov)
    }

    // Per-dataset h5ad sources tagged with the dataset name: tuple(dataset, h5ad).
    ds_meta = channel.fromList(
        active_configs.collect { name, dcfg -> tuple(name, file(dcfg.sc_data)) }
    )

    // Per-dataset VCFs tagged with the dataset name: tuple(dataset, chr, vcf).
    vcf_files = channel.empty()
    active_configs.each { name, dcfg ->
        vcf_files = vcf_files.mix(
            channel.fromFilePairs(dcfg.vcf_glob, size: 1, flat: true)
                .map { chr, vcf -> tuple(name, chr, vcf) }
        )
    }

    ids = EXTRACT_INDIV_IDS(ds_meta)
    cluster_sizes = COMPUTE_CLUSTER_SIZES(ds_meta)
    filt_vcf_files = FILTER_VCF(vcf_files.combine(ids, by: 0))
    bed_files = CONVERT_VCF_TO_BED(filt_vcf_files)
        .map { dataset, chr, bed, bim, fam -> tuple(dataset, chr, bed) }
    pruned_snps = PRUNE_SNPS(bed_files)
    all_bed = CONCAT_BED_FILES(
        bed_files.map { dataset, chr, bed -> tuple(dataset, bed) }.groupTuple()
    ).map { dataset, bed, bim, fam -> tuple(dataset, bed) }
    geno_pcs = COMPUTE_GENOTYPE_PCS(all_bed)
    grm = CREATE_GRM(all_bed)
    anno_bed = CREATE_ANNOT_BED(params.annot_bed_url)

    // Per-dataset cell_infos; info.dataset is set on every info map.
    cell_infos = channel.empty()
    active_configs.each { name, dcfg ->
        def ds_infos = channel.fromList(dcfg.cell_types)
            .flatMap { ct ->
                def info = [dataset: name, cell_type: ct.toString(), cell_frac: 1.0d, indiv_frac: 1.0d, n_cells_target: -1, count_frac: 1.0d]
                def combos = []
                if (info.cell_type in dcfg.subsample_cell_types) {
                    def cell_fracs = [0.005d, 0.01d, 0.05d, 0.1d, 0.25d, 0.5d, 1.0d]
                    combos += cell_fracs.collect { f -> info + [cell_frac: f, indiv_frac: 1.0d, count_frac: 1.0d] }
                    if (info.cell_type in dcfg.indiv_frac_cell_types) {
                        def indiv_fracs = [0.05d, 0.1d, 0.25d, 0.5d, 1.0d]
                        combos += indiv_fracs.collect { f -> info + [cell_frac: 1.0d, indiv_frac: f, count_frac: 1.0d] }
                        // Joint cell × individual grid for CD4_NC power heatmap.
                        combos += joint_fracs.collectMany { cf ->
                            joint_fracs.collect { inf ->
                                info + [cell_frac: cf, indiv_frac: inf, n_cells_target: -1, count_frac: 1.0d, joint_both_frac: true]
                            }
                        }
                    }
                    def n_cells_targets = [CD4_NC: 300, B_IN: 100]
                    if (n_cells_targets.containsKey(info.cell_type)) {
                        combos += [[
                            dataset: name,
                            cell_type: info.cell_type,
                            cell_frac: 1.0d,
                            indiv_frac: 1.0d,
                            n_cells_target: n_cells_targets[info.cell_type],
                            count_frac: 1.0d
                        ]]
                    }
                    // Independent count-depth series on full cells / individuals.
                    def count_fracs = [0.1d, 0.25d, 0.5d, 1.0d]
                    combos += count_fracs.collect { f ->
                        info + [cell_frac: 1.0d, indiv_frac: 1.0d, n_cells_target: -1, count_frac: f]
                    }
                    combos = combos.unique { x -> "${x.cell_frac}|${x.indiv_frac}|${x.n_cells_target}|${x.count_frac}" }
                    return combos
                }
                return [info]
            }

        if (!dcfg.combined_types.isEmpty()) {
            ds_infos = ds_infos.mix(
                channel.fromList(
                    dcfg.combined_types.collect { cname, members ->
                        [dataset: name, cell_type: cname, cell_frac: 1.0d, indiv_frac: 1.0d, n_cells_target: -1, count_frac: 1.0d]
                    }
                )
            )
        }
        cell_infos = cell_infos.mix(ds_infos)
    }

    cell_infos = cell_infos.map { info ->
        def has_n_cells_filter = (info.n_cells_target >= 0)
        def is_full_data = (info.cell_frac == 1.0d && info.indiv_frac == 1.0d && !has_n_cells_filter && info.count_frac == 1.0d)
        info + [is_full_data: is_full_data]
    }

    // TEMP: disable joint and other downsampling — keep full data only.
    cell_infos = cell_infos.filter { it.is_full_data }

    // Pair each info with its dataset's h5ad.
    sc_input = cell_infos
        .map { info -> tuple(info.dataset, info) }
        .combine(ds_meta, by: 0)
        .map { dataset, info, scf -> tuple(info, scf) }

    gene_properties = COMPUTE_GENE_PROPERTIES(sc_input, anno_bed)
    gene_properties_plot = gene_properties.filter { info, rest -> info.is_full_data }
    gene_properties_file = Channel
        .of("dataset\tcell_type\tproperties_file")
        .concat(gene_properties_plot
            .map { info, properties_file->
                "${info.dataset}\t${info.cell_type}\t${properties_file}"
            }
        )
        .collectFile(name: 'gene_properties_file', newLine: true, sort: false)

    sc_input = sc_input.combine(gene_properties, by: 0)
    pb_counts = COMPUTE_PB_COUNTS(sc_input)
    pb_logcounts = COMPUTE_PB_LOGCOUNTS(sc_input)
    pb_input = pb_counts.mix(pb_logcounts)
    pheno = ANNOTATE_PHENO(pb_input, anno_bed)
    sc_counts = COMPUTE_SC_COUNTS(sc_input)
    sc_sct_counts = COMPUTE_SC_SCT_COUNTS(
        sc_counts
            .filter { _sc_type, info, _pheno ->
                info.cell_type == "B_all" && info.is_full_data
            }
            .map { _sc_type, info, pheno -> tuple(info.dataset, info, pheno) }
            .combine(ds_meta, by: 0)
            .map { _dataset, info, pheno, h5ad -> tuple(info, pheno, h5ad) }
    )
    sc_logcounts = COMPUTE_SC_LOGCOUNTS(
        sc_input.filter { it[0].is_full_data
                       || (it[0].indiv_frac == 1.0d && it[0].count_frac == 1.0d)
                       || (it[0].cell_type == "CD4_NC" && it[0].cell_frac == 1.0d && it[0].count_frac == 1.0d)
                       || is_joint_both_frac(it[0])
                       || (it[0].n_cells_target >= 0 && it[0].count_frac == 1.0d)
                       || (it[0].cell_frac == 1.0d
                           && it[0].indiv_frac == 1.0d
                           && it[0].n_cells_target < 0
                           && it[0].count_frac != 1.0d) }
    )
    sc_pheno = sc_counts.mix(sc_logcounts).mix(sc_sct_counts)

    pb_expr_covs = COMPUTE_PB_EXPR_COVS(sc_input)
    sc_expr_covs = COMPUTE_SC_EXPR_COVS(sc_input)
    starcat_ref = DOWNLOAD_STARCAT_REF(params.starcat_ref_url)
    starcat_covs = COMPUTE_STARCAT_COVS(
        sc_input.filter {it[0].cell_type in ["T_all"]}
            .combine(starcat_ref)
    )

    slingshot_input = cell_infos
        .filter { it.cell_type == "B_all" }
        .map { info -> tuple(info.dataset, info) }
        .combine(ds_meta, by: 0)
        .map { dataset, info, scf -> tuple(info, scf) }
    slingshot_adata = PREPARE_SLINGSHOT_ADATA(slingshot_input)
    slingshot_out = RUN_SLINGSHOT(slingshot_adata)

    seacells_input = cell_infos
        .filter { it.is_full_data }
        .filter { !(it.cell_type in ["Plasma", "CD4_SOX4"]) }
        .map { info -> tuple(info.dataset, info) }
        .combine(ds_meta, by: 0)
        .map { dataset, info, scf -> tuple(info, scf) }
    seacells_adata = PREPARE_SEACELLS_ADATA(seacells_input)
    // TEMP: skip SEACells on T_all (too many cells).
    seacells_out = COMPUTE_SEACELLS(seacells_adata.filter { it[0].cell_type != "T_all" })
    seacells_groups = seacells_out.groups

    sc_preprocess_input = sc_input
    sc_input = sc_counts
        .map { sc_type, info, sc_counts_file -> tuple(info, sc_counts_file) }
        .combine(pb_expr_covs, by: 0)
        .map { info, sc_counts_file, pb_covs_file -> tuple(info.dataset, info, sc_counts_file, pb_covs_file) }
        .combine(geno_pcs, by: 0)
        .map { dataset, info, sc_counts_file, pb_covs_file, geno_pcs_file ->
            tuple(info, sc_counts_file, pb_covs_file, geno_pcs_file)
        }

    // Run SAIGE-QTL. Skip T_all (too many cells for collation).
    saige_subset_bed = SAIGE_SUBSET_BED(all_bed)
    saigeqtl_input = sc_input
        .filter { it[0].is_full_data && it[0].cell_type != "T_all" }
        .map { info, sc_counts_file, pb_covs_file, geno_pcs_file ->
            tuple(info.dataset, info, sc_counts_file, pb_covs_file, geno_pcs_file)
        }
        .combine(bed_files, by: 0)
        .map { dataset, info, sc_counts_file, pb_covs_file, geno_pcs_file, chr, bed_file ->
            tuple(info + [chr: chr], sc_counts_file, pb_covs_file, geno_pcs_file, bed_file)
        }
        .combine(anno_bed)
        .map { info, sc_counts_file, pb_covs_file, geno_pcs_file, bed_file, anno ->
            tuple(info, sc_counts_file, pb_covs_file, geno_pcs_file, bed_file, anno)
        }
        // TEMP: disable SAIGE-QTL.
        .filter { false }
        | COLLATE_SAIGEQTL_INPUT

    gene_lists = EXTRACT_GENES(saigeqtl_input.map({[it[0], it[1]]}), 50)
        .flatMap { info, gene_lists ->
            def lists = (gene_lists instanceof List) ? gene_lists : [gene_lists]
            lists.collect { glist -> tuple(info, glist) }
        }

    saigeqtl_input = saigeqtl_input 
        .combine(gene_lists, by: 0)
        .filter { it[0].is_full_data }
        .map { it -> [it[0].dataset] + it }
        .combine(saige_subset_bed, by: 0)
        .map { dataset, info, input, bed, annot, gene_file, subset_bed ->
            tuple(info, input, bed, annot, gene_file, subset_bed)
        }
        .filter { it[0].cell_type == "Plasma" }

    saigeqtl_output = RUN_SAIGEQTL(saigeqtl_input)
    saigeqtl_output = COMPUTE_SAIGEQTL_POWER(saigeqtl_output)
    // TEMP: skip SAIGE-QTL concat while SAIGE-QTL is disabled (empty collect never emits).
    saigeqtl_file = Channel.empty()
    // saigeqtl_file = CONCAT_SAIGEQTL_TSVS(saigeqtl_output.map{it[1]}.collect())

    // Run pseudobulk quasar.
    pb_covs = COMBINE_PB_COVS(
        pb_expr_covs
            .map { info, covs -> tuple(info.dataset, info, covs) }
            .combine(geno_pcs, by: 0)
            .map { dataset, info, covs, pcs -> tuple(info, covs, pcs) }
    )

    pc_gwas_input = COMPUTE_SC_PC_PHENO(
        sc_preprocess_input
            .filter { it[0].is_full_data && !(it[0].cell_type in ["B_all", "T_all"]) }
            .map { info, scf, _gene_props -> tuple(info, scf) }
    )
        .map { info, pheno -> tuple([info.dataset, info.cell_type], info, pheno) }
        .combine(
            pb_covs
                .filter { it[0].is_full_data && !(it[0].cell_type in ["B_all", "T_all"]) }
                .map { info, covs -> tuple([info.dataset, info.cell_type], covs) },
            by: 0
        )
        .map { key, info, pheno, covs ->
            tuple(info.dataset, info + [model: 'lm', gwas_label: 'pc'], pheno, covs)
        }
        .combine(all_bed, by: 0)
        .map { dataset, info, pheno, covs, bed -> tuple(info, pheno, covs, bed) }

    pc_gwas = RUN_PC_GWAS(pc_gwas_input)

    int_cov = Channel.of("none", "age", "sex")
    models_for_pb_type = ["counts": ["nb_glm", "nb_glmm", "p_glmm"], "logcounts": ["lm", "lmm"]]
    pb_quasar_input = pheno
        .combine(pb_covs, by: 0)
        .map { info, pb_type_val, pheno_bed, covs -> tuple(info.dataset, info, pb_type_val, pheno_bed, covs) }
        .combine(bed_files, by: 0)
        .combine(grm, by: 0)
        .combine(int_cov)
        .flatMap { dataset, info, pb_type_val, pheno_bed, covs, chr, bed_file, grm_path, int_cov_val ->
            models_for_pb_type[pb_type_val].collect { model_val ->
                def dict = info + [chr: chr, model: model_val, int_cov: int_cov_val, pb_type: pb_type_val]
                tuple(dict, pheno_bed, covs, bed_file, grm_path)
            }
        }
        .filter({
            it[0].is_full_data ||
            (it[0].model in ["lm", "nb_glm"] && it[0].indiv_frac == 1.0d && it[0].count_frac == 1.0d) ||
            (it[0].model == "lm" && it[0].cell_type == "CD4_NC" && it[0].cell_frac == 1.0d && it[0].count_frac == 1.0d) ||
            (it[0].model == "lm" && is_joint_both_frac(it[0])) ||
            (it[0].model in ["lm", "nb_glm"] && it[0].n_cells_target >= 0 && it[0].count_frac == 1.0d) ||
            (it[0].model in ["lm", "nb_glm"]
                && it[0].cell_frac == 1.0d
                && it[0].indiv_frac == 1.0d
                && it[0].n_cells_target < 0
                && it[0].count_frac != 1.0d)
        })
        .filter({ !(it[0].cell_type in ["B_all", "T_all"]) })
        // TEMP: all PB models, B_IN / Plasma Gxsex only (nb_glmm: Plasma only).
        .filter({
            it[0].cell_type in ["B_IN", "Plasma"] &&
            it[0].int_cov == "sex" &&
            it[0].is_full_data &&
            !(it[0].model == "nb_glmm" && it[0].cell_type == "B_IN")
        })

    pb_quasar = RUN_QUASAR_PB(pb_quasar_input)
    pb_quasar = Utils.attachGeneProperties(pb_quasar, gene_properties)

    pb_quasar = COMPUTE_QUASAR_POWER_PB(pb_quasar)
    pb_quasar = COMPUTE_CONVERGENCE_PB(pb_quasar)

    // Run single-cell quasar.
    starcat_int_cov_by_cell_type = [
        "T_all": "starcat_all",
    ]
    starcat_cell_types = starcat_int_cov_by_cell_type.keySet() as List
    starcat_int_covs = starcat_int_cov_by_cell_type.values() as List

    sc_expr_pseudo_input = sc_expr_covs
        .filter { it[0].cell_type == "B_all" }
        .map { info, expr -> tuple(info.dataset, info, expr) }
        .combine(slingshot_out.map { sinfo, pt, plot1, plot2 -> tuple(sinfo.dataset, pt) }, by: 0)
        .map { dataset, info, expr, pt -> tuple(info, expr, pt) }
    sc_expr_starcat_input = sc_expr_covs
        .filter { starcat_cell_types.contains(it[0].cell_type) }
        .combine(starcat_covs, by: 0)
        .map { info, expr, starcat_cov, starcat_plot -> tuple(info, expr, starcat_cov) }
    sc_expr_aug = JOIN_SC_INT_COVS(sc_expr_pseudo_input.mix(sc_expr_starcat_input))
    sc_expr_covs_aug = sc_expr_covs
        .filter { !(it[0].cell_type == "B_all" || starcat_cell_types.contains(it[0].cell_type)) }
        .mix(sc_expr_aug)

    both_expr_covs = pb_expr_covs.combine(sc_expr_covs_aug, by: 0)
    sc_covs_combined = COMBINE_SC_COVS(
        both_expr_covs
            .map { info, pb_c, sc_c -> tuple(info.dataset, info, pb_c, sc_c) }
            .combine(geno_pcs, by: 0)
            .map { dataset, info, pb_c, sc_c, pcs -> tuple(info, pb_c, sc_c, pcs) }
    )

    starcat_cov_names = [
        "starcat_Cytotoxic",
        "starcat_TEMRA",
        "starcat_CD4_CM",
        "starcat_CD8_EM",
        "starcat_CD4_Naive",
    ]
    base_cov_specs = Channel.of("bulk_pca", "sc_pca", "bulk_pca+pct_mito", "bulk_pca+cell_cycle")
    int_cov_specs = Channel.fromList([
        ["B_all", "bulk_pca+pseudotime"],
        ["T_all", "bulk_pca+starcat_all"],
        ["T_all", "bulk_pca+starcat_CD4_Naive"],
    ])
    int_cov_for_cov_spec = [
        "bulk_pca+pseudotime": "pseudotime",
        "bulk_pca+starcat_all": "starcat_all",
        "bulk_pca+starcat_CD4_Naive": "starcat_CD4_Naive",
    ]
    def interaction_cov_for_int_cov = { int_cov ->
        int_cov == "starcat_all" ? starcat_cov_names.join(",") : int_cov
    }
    base_cov_pairs = sc_covs_combined.combine(base_cov_specs)
    int_cov_pairs = sc_covs_combined
        .map { info, cov -> tuple(info.cell_type, info, cov) }
        .combine(int_cov_specs, by: 0)
        .map { cell_type, info, cov, cov_spec -> tuple(info, cov, cov_spec) }
    sc_covs = FILTER_COVS(base_cov_pairs.mix(int_cov_pairs))

    // TEMP: CASTIE on T_all counts with starcat_CD4_Naive only.
    castie_pheno = sc_counts
        .filter { it[1].is_full_data && it[1].cell_type == "T_all" }
        .map { data_type, info, pheno_file -> tuple(info, data_type, pheno_file) }
    castie_covs = sc_covs.filter { info, cov_spec, _covs ->
        info.is_full_data &&
            info.cell_type == "T_all" &&
            cov_spec == "bulk_pca+starcat_CD4_Naive"
    }
    castie_input = castie_pheno
        .combine(castie_covs, by: 0)
        .map { info, data_type, pheno_file, cov_spec, covs ->
            tuple(info.dataset, info, data_type, pheno_file, covs, cov_spec)
        }
        .combine(bed_files, by: 0)
        .map { dataset, info, data_type, pheno_file, covs, cov_spec, chr, bed_file ->
            tuple(
                info + [
                    chr: chr,
                    data_type: data_type,
                    cov_spec: cov_spec,
                    int_cov: int_cov_for_cov_spec.getOrDefault(cov_spec, "none"),
                ],
                pheno_file, covs, bed_file
            )
        }
        .combine(anno_bed)
        .map { info, pheno_file, covs, bed_file, anno ->
            tuple(info, pheno_file, covs, bed_file, anno)
        }
        // TEMP: CASTIE run disabled.
        // | COLLATE_CASTIE_INPUT

    // castie_gene_lists = EXTRACT_CASTIE_GENES(castie_input.map({[it[0], it[1]]}), 10)
    //     .flatMap { info, gene_lists ->
    //         def lists = (gene_lists instanceof List) ? gene_lists : [gene_lists]
    //         lists.collect { glist -> tuple(info, glist) }
    //     }

    // castie_run_input = castie_input
    //     .combine(castie_gene_lists, by: 0)
    //     .map { it -> [it[0].dataset] + it }
    //     .combine(saige_subset_bed, by: 0)
    //     .map { dataset, info, input, bed, annot, gene_file, subset_bed ->
    //         tuple(info, input, bed, annot, gene_file, subset_bed)
    //     }

    // castie_output = RUN_CASTIE(castie_run_input)
    // castie_file = CONCAT_CASTIE_TSVS(castie_output.map{it[1]}.collect())

    pc_sc_gwas_input = COMPUTE_SC_PC_SC_PHENO(
        sc_preprocess_input
            .filter { it[0].is_full_data && !(it[0].cell_type in ["B_all", "T_all"]) }
            .map { info, scf, _gene_props -> tuple(info, scf) }
    )
        .map { info, pheno -> tuple([info.dataset, info.cell_type], info, pheno) }
        .combine(
            sc_covs
                .filter { it[0].is_full_data && it[1] == "bulk_pca" && !(it[0].cell_type in ["B_all", "T_all"]) }
                .map { info, cov_spec, covs -> tuple([info.dataset, info.cell_type], covs) },
            by: 0
        )
        .map { key, info, pheno, covs ->
            tuple(info.dataset, info + [model: 'lmm_sc', gwas_label: 'pc'], pheno, covs)
        }
        .combine(all_bed, by: 0)
        .map { dataset, info, pheno, covs, bed -> tuple(info, pheno, covs, bed) }

    pc_sc_gwas = RUN_PC_SC_GWAS(pc_sc_gwas_input)

    models_for_data_type = [
        "counts": ["p_glmm_sc"],
        "log_counts": ["lmm_sc"],
        "sct_counts": ["p_glmm_sc"],
    ]
    sc_quasar_input = sc_pheno
        .map { data_type, info, sc_pheno_file -> tuple(info, data_type, sc_pheno_file) }
        .combine(sc_covs, by: 0)
        .combine(anno_bed)
        .map { info, data_type, sc_pheno_file, cov_spec, covs, anno ->
            tuple(info.dataset, info, data_type, sc_pheno_file, cov_spec, covs, anno)
        }
        .combine(bed_files, by: 0)
        .flatMap { dataset, info, data_type, sc_pheno_file, cov_spec, covs, anno, chr, bed_file ->
            models_for_data_type[data_type].collect { model_val ->
                tuple(
                    info + [
                        chr: chr,
                        cov_spec: cov_spec,
                        data_type: data_type,
                        model: model_val,
                        int_cov: int_cov_for_cov_spec.getOrDefault(cov_spec, "none"),
                        interaction_cov: interaction_cov_for_int_cov(
                            int_cov_for_cov_spec.getOrDefault(cov_spec, "none")
                        ),
                    ],
                    sc_pheno_file, covs, anno, bed_file
                )
            }
        }
        .filter { !(it[0].cov_spec in ["sc_pca", "bulk_pca+pct_mito", "bulk_pca+cell_cycle"])
                || (it[0].cell_type in ["Plasma", "B_IN", "CD4_NC"]) }
        .filter { it[0].cov_spec == "bulk_pca" || it[0].model == "p_glmm_sc" || it[0].int_cov != "none" }
        .filter { it[0].is_full_data
               || (it[0].model in ["p_glmm_sc", "lmm_sc"]
                   && it[0].cov_spec == "bulk_pca"
                   && it[0].indiv_frac == 1.0d
                   && it[0].count_frac == 1.0d)
               || (it[0].model in ["p_glmm_sc", "lmm_sc"]
                   && it[0].cov_spec == "bulk_pca"
                   && it[0].cell_type == "CD4_NC"
                   && it[0].cell_frac == 1.0d
                   && it[0].count_frac == 1.0d)
               || (it[0].model == "p_glmm_sc"
                   && it[0].cov_spec == "bulk_pca"
                   && is_joint_both_frac(it[0]))
               || (it[0].model in ["p_glmm_sc", "lmm_sc"]
                   && it[0].cov_spec == "bulk_pca"
                   && it[0].n_cells_target >= 0
                   && it[0].count_frac == 1.0d)
               || (it[0].model in ["p_glmm_sc", "lmm_sc"]
                   && it[0].cov_spec == "bulk_pca"
                   && it[0].cell_frac == 1.0d
                   && it[0].indiv_frac == 1.0d
                   && it[0].n_cells_target < 0
                   && it[0].count_frac != 1.0d) }

    sc_cell_group_k = Channel.of("none", "seacells", 5)
    sc_input_by_mode = sc_quasar_input
        .combine(sc_cell_group_k)
        .map { info, sc_pheno, covs, anno, plink, k -> tuple(info + [k: k], sc_pheno, covs, anno, plink) }
        .branch {
            none:     it[0].k == "none"
            seacells: it[0].k == "seacells"
            binned:   it[0].k != "none" && it[0].k != "seacells"
            drop:     true
        }

    sc_none_full = sc_input_by_mode.none
        .map { info, sc_pheno, covs, anno, plink -> tuple(info, sc_pheno, covs, anno, plink, covs) }

    sc_seacells_full = sc_input_by_mode.seacells
        .map { info, sc_pheno, covs, anno, plink ->
            def base = [dataset: info.dataset, cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac]
            tuple(base, info, sc_pheno, covs, anno, plink)
        }
        .combine(
            seacells_groups.map { ginfo, cg ->
                tuple([dataset: ginfo.dataset, cell_type: ginfo.cell_type, cell_frac: ginfo.cell_frac, indiv_frac: ginfo.indiv_frac], cg)
            },
            by: 0
        )
        .map { base, info, sc_pheno, covs, anno, plink, cg ->
            tuple(info, sc_pheno, covs, anno, plink, cg)
        }

    def allow_grouped_int_cov = { int_cov ->
        int_cov != "none" && (!is_starcat_int_cov(int_cov) || int_cov == "starcat_CD4_Naive")
    }
    sc_binned_groups = COMPUTE_CELL_GROUPS(
        sc_input_by_mode.binned
            .filter { allow_grouped_int_cov(it[0].int_cov) }
            .map { info, sc_pheno, covs, anno, plink -> tuple(info, covs) }
    )
    sc_binned_full = sc_input_by_mode.binned
        .filter { allow_grouped_int_cov(it[0].int_cov) }
        .combine(sc_binned_groups, by: 0)
        .map { info, sc_pheno, covs, anno, plink, cg ->
            tuple(info, sc_pheno, covs, anno, plink, cg)
        }

    sc_quasar_input_full = sc_none_full
        .mix(sc_seacells_full)
        .mix(sc_binned_full)
        .filter({ is_int_sc_cov(it[0].int_cov) })
        .filter({it[0].k in ["none", 5]})
        .filter({ !is_starcat_int_cov(it[0].int_cov) || it[0].k == "none" || it[0].int_cov == "starcat_CD4_Naive" })
        .filter { it[0].model != "lmm_sc"
               || (it[0].dataset == "onek1k"
                   && it[0].k == "none") }
        .filter { it[0].cell_type in ["B_all", "T_all"] }
        .filter { it[0].data_type != "sct_counts"
               || (it[0].cell_type == "B_all"
                   && it[0].int_cov in ["none", "pseudotime"]
                   && it[0].k == "none") }
        .filter { it[0].is_full_data
               || (it[0].model == "p_glmm_sc"
                   && it[0].k == 5
                   && it[0].int_cov != "none"
                   && it[0].cell_type in ["B_all"])
               || (it[0].model in ["p_glmm_sc", "lmm_sc"]
                   && it[0].k == "none"
                   && it[0].int_cov == "none"
                   && it[0].indiv_frac == 1.0d
                   && it[0].count_frac == 1.0d)
               || (it[0].model in ["p_glmm_sc", "lmm_sc"]
                   && it[0].k == "none"
                   && it[0].int_cov == "none"
                   && it[0].cell_type == "CD4_NC"
                   && it[0].cell_frac == 1.0d
                   && it[0].count_frac == 1.0d)
               || (it[0].model == "p_glmm_sc"
                   && it[0].k == "none"
                   && it[0].int_cov == "none"
                   && is_joint_both_frac(it[0]))
               || (it[0].model in ["p_glmm_sc", "lmm_sc"]
                   && it[0].k == "none"
                   && it[0].int_cov == "none"
                   && it[0].n_cells_target >= 0
                   && it[0].count_frac == 1.0d)
               || (it[0].model in ["p_glmm_sc", "lmm_sc"]
                   && it[0].k == "none"
                   && it[0].int_cov == "none"
                   && it[0].cell_frac == 1.0d
                   && it[0].indiv_frac == 1.0d
                   && it[0].n_cells_target < 0
                   && it[0].count_frac != 1.0d) }
        // TEMP: B_all pseudotime + T_all starcat interaction runs only.
        .filter {
            (it[0].model == "p_glmm_sc"
                && it[0].cell_type == "B_all"
                && it[0].data_type == "sct_counts"
                && it[0].k == "none"
                && it[0].int_cov == "pseudotime")
            || (it[0].model == "p_glmm_sc"
                && it[0].cell_type == "T_all"
                && it[0].data_type == "counts"
                && it[0].k == "none"
                && it[0].int_cov in ["starcat_all", "starcat_CD4_Naive"])
        }

    sc_quasar = RUN_QUASAR_SC(sc_quasar_input_full)
    sc_quasar = Utils.attachGeneProperties(sc_quasar, gene_properties)

    sc_quasar = COMPUTE_QUASAR_POWER_SC(sc_quasar)
    sc_quasar = COMPUTE_CONVERGENCE_SC(sc_quasar)

    // Offset-precision experiment. Everything except the offset supplied to
    // p_glmm_sc is held fixed, so any power difference is attributable to how
    // precisely sequencing depth is described.
    sc_offset_base = sc_counts
        .map { sc_type, info, sc_pheno_file -> tuple(info, sc_pheno_file) }
        .filter { it[0].cell_type == "B_IN" && it[0].is_full_data }

    sc_offset_files = COMPUTE_SC_OFFSETS(sc_offset_base)
        .flatMap { info, percell, donorflat, constant, summary ->
            [tuple(info, "percell", percell),
             tuple(info, "donorflat", donorflat),
             tuple(info, "constant", constant)]
        }

    sc_offset_quasar_input = sc_offset_base
        .combine(sc_covs.filter { it[1] == "bulk_pca" }, by: 0)
        .combine(sc_offset_files, by: 0)
        .combine(anno_bed)
        .map { info, sc_pheno_file, cov_spec, covs, offset_spec, offset_file, anno ->
            tuple(info.dataset, info, sc_pheno_file, cov_spec, covs, offset_spec, offset_file, anno)
        }
        .combine(bed_files, by: 0)
        .map { dataset, info, sc_pheno_file, cov_spec, covs, offset_spec, offset_file, anno, chr, bed_file ->
            def dict = info + [
                chr: chr,
                cov_spec: cov_spec,
                offset_spec: offset_spec,
                data_type: "counts",
                model: "p_glmm_sc",
                int_cov: "none",
                k: "none",
            ]
            tuple(dict, sc_pheno_file, covs, offset_file, anno, bed_file)
        }

    sc_offset_quasar = RUN_QUASAR_SC_OFFSET(sc_offset_quasar_input)
    sc_offset_quasar = Utils.attachGeneProperties(sc_offset_quasar, gene_properties)
    sc_offset_quasar = COMPUTE_QUASAR_POWER_SC_OFFSET(sc_offset_quasar)

    // Run csqQTL analysis.
    seacells_csaqtl_groups = seacells_groups
        .map { info, cg ->
            tuple([dataset: info.dataset, cell_type: info.cell_type, cell_frac: info.cell_frac,
                   indiv_frac: info.indiv_frac, grouping: "seacells"], cg)
        }

    int_csaqtl_groups_input = sc_covs
        .filter { info, cov_spec, cov ->
            int_cov_for_cov_spec.containsKey(cov_spec) &&
                !is_starcat_int_cov(int_cov_for_cov_spec[cov_spec])
        }
        .map { info, cov_spec, cov ->
            tuple(info + [int_cov: int_cov_for_cov_spec[cov_spec], k: 5], cov)
        }
    int_csaqtl_groups = COMPUTE_CSAQTL_GROUPS(int_csaqtl_groups_input)
        .map { info, cg ->
            tuple([dataset: info.dataset, cell_type: info.cell_type, cell_frac: info.cell_frac,
                   indiv_frac: info.indiv_frac, grouping: info.int_cov], cg)
        }
    csaqtl_cell_groups = seacells_csaqtl_groups.mix(int_csaqtl_groups)

    csaqtl_pheno = COMPUTE_CSAQTL_COUNTS(
        csaqtl_cell_groups
            .map { info, cg -> tuple([info.dataset, info.cell_type], info, cg) }
            .combine(seacells_adata.map { ainfo, h5ad -> tuple([ainfo.dataset, ainfo.cell_type], h5ad) }, by: 0)
            .map { key, info, cg, h5ad -> tuple(info, cg, h5ad) }
    )

    csaqtl_quasar_input = csaqtl_pheno
        .map { info, pheno -> tuple([info.dataset, info.cell_type], info, pheno) }
        .combine(
            pb_covs.filter { it[0].is_full_data }.map { pinfo, covs -> tuple([pinfo.dataset, pinfo.cell_type], covs) },
            by: 0
        )
        .map { key, info, pheno, covs -> tuple(info.dataset, info + [model: 'nb_glm', gwas_label: 'csaqtl'], pheno, covs) }
        .combine(all_bed, by: 0)
        .map { dataset, info, pheno, covs, bed -> tuple(info, pheno, covs, bed) }

    csaqtl_quasar = RUN_CSAQTL_GWAS(csaqtl_quasar_input)

    // Plotting output (int figures currently disabled; extracts need gene/SNP list TSVs).
    // sc_logcounts = EXTRACT_SC_LOGCOUNTS(sc_preprocess_input.filter { it[0].is_full_data })
    // genotype_dosages = EXTRACT_GENOTYPES(all_bed)
    //
    // sc_logcounts_for_figures = sc_logcounts
    //     .map { info, logcounts -> tuple([info.dataset, info.cell_type], logcounts) }
    // sc_int_figures_input = sc_expr_aug
    //     .filter { it[0].cell_type == "B_all" }
    //     .map { info, covs -> tuple([info.dataset, info.cell_type], info, covs) }
    //     .combine(sc_logcounts_for_figures, by: 0)
    //     .map { key, info, covs, logcounts -> tuple(info.dataset, info, covs, logcounts) }
    //     .combine(genotype_dosages, by: 0)
    //     .map { dataset, info, covs, logcounts, geno ->
    //         tuple(info, covs, logcounts, geno)
    //     }

    pb_quasar_file = Channel
        .of("dataset\tmodel\tcell_type\tchr\tcell_frac\tindiv_frac\tn_cells_target\tcount_frac\tint_cov\tpb_type\tregion_file\tvariant_file\ttime_file\tpower_file\tconv_file\tgene_prop_file")
        .concat(pb_quasar
            .map { info, region_file, variant_file, time_file, gene_prop_file, power_file, conv_file  ->
                "${info.dataset}\t${info.model}\t${info.cell_type}\t${info.chr}\t${info.cell_frac}\t${info.indiv_frac}\t${info.n_cells_target}\t${info.count_frac}\t${info.int_cov}\t${info.pb_type}\t${region_file}\t${variant_file}\t${time_file}\t${power_file}\t${conv_file}\t${gene_prop_file}"
            }
        )
        .collectFile(name: 'pb_quasar_file', newLine: true, sort: false)
    
    sc_quasar_file = Channel
        .of("dataset\tmodel\tdata_type\tcell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tn_cells_target\tcount_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tpower_file\tconv_file\tgene_prop_file")
        .concat(sc_quasar
            .map { info, region_file, variant_file, time_file, gene_prop_file, power_file, conv_file ->
                "${info.dataset}\t${info.model}\t${info.data_type}\t${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.n_cells_target}\t${info.count_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${power_file}\t${conv_file}\t${gene_prop_file}"
            }
        )
        .collectFile(name: 'sc_quasar_file', newLine: true, sort: false)

    sc_offset_quasar_file = Channel
        .of("dataset\tmodel\tcell_type\tchr\tcov_spec\toffset_spec\tregion_file\tvariant_file\ttime_file\tpower_file\tgene_prop_file")
        .concat(sc_offset_quasar
            .map { info, region_file, variant_file, time_file, gene_prop_file, power_file ->
                "${info.dataset}\t${info.model}\t${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.offset_spec}\t${region_file}\t${variant_file}\t${time_file}\t${power_file}\t${gene_prop_file}"
            }
        )
        .collectFile(name: 'sc_offset_quasar_file', newLine: true, sort: false)

    seacells_file = Channel
        .of("dataset\tcell_type\tgroups_file\tinfo_file")
        .concat(seacells_groups
            .combine(seacells_out.sizes, by: 0)
            .map { info, cg, sizes ->
                "${info.dataset}\t${info.cell_type}\t${cg}\t${sizes}"
            }
        )
        .collectFile(name: 'seacells_file', newLine: true, sort: false)

    // Clumping analysis.
    // clumped_pb_quasar = CLUMP_VARIANTS_PB(
    //     pb_quasar.filter{it[0].int_cov == "none"}
    //         .map { it -> [it[0].dataset] + it }
    //         .combine(all_bed, by: 0)
    //         .map { it -> it[1..-1] }
    // )
    // clumped_pb_quasar_file = Channel
    //     .of("dataset\tmodel\tcell_type\tchr\tcell_frac\tindiv_frac\tint_cov\tpb_type\tregion_file\tvariant_file\ttime_file\tpower_file\tconv_file\tgene_prop_file\tclumped_file")
    //     .concat(clumped_pb_quasar
    //         .map { info, region_file, variant_file, time_file, gene_prop_file, power_file, conv_file, clumped_file ->
    //             "${info.dataset}\t${info.model}\t${info.cell_type}\t${info.chr}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.pb_type}\t${region_file}\t${variant_file}\t${time_file}\t${power_file}\t${conv_file}\t${gene_prop_file}\t${clumped_file}"
    //         }
    //     )
    //     .collectFile(name: 'clumped_pb_quasar_file', newLine: true, sort: false)

    // clumped_sc_quasar = CLUMP_VARIANTS_SC(
    //     sc_quasar.filter{it[0].int_cov == "none"}
    //         .map { it -> [it[0].dataset] + it }
    //         .combine(all_bed, by: 0)
    //         .map { it -> it[1..-1] }
    // )
    // clumped_sc_quasar_file = Channel
    //     .of("dataset\tcell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tpower_file\tconv_file\tgene_prop_file\tclumped_file")
    //     .concat(clumped_sc_quasar
    //         .map { info, region_file, variant_file, time_file, gene_prop_file, power_file, conv_file, clumped_file ->
    //             "${info.dataset}\t${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${power_file}\t${conv_file}\t${gene_prop_file}\t${clumped_file}"
    //         }
    //     )
    //     .collectFile(name: 'clumped_sc_quasar_file', newLine: true, sort: false)

    // GWAS analysis.
    // Disabled for now.
    // sc_quasar_gwas_input = sc_quasar_input
    //     .filter({it[0].model == "lmm_sc"})
    //     .filter({it[0].cov_spec == "bulk_pca" && it[0].is_full_data})
    //     .filter({it[0].cell_type == "B_IN"})
    //     .combine(Channel.from(1..22))
    //     .map { info, sc_pheno, covs, anno, bed, pheno_chr ->
    //       tuple(info + [pheno_chr: pheno_chr as int], sc_pheno, covs, anno, bed)
    //     }
    //     .filter({it[0].pheno_chr== 1.0d})
    // sc_quasar_gwas = RUN_QUASAR_SC_GWAS(sc_quasar_gwas_input)
    sc_quasar_gwas = Channel.empty()

    pb_quasar_gwas_input = pb_quasar_input
        .filter({it[0].is_full_data && it[0].int_cov == "none" && it[0].model == "nb_glm"})
        .filter({it[0].cell_type == "B_IN"})
        .combine(Channel.from(1..22))
        .map { info, pheno, covs, bed, grm, pheno_chr ->
          tuple(info + [pheno_chr: pheno_chr as int], pheno, covs, bed)
        }
        .filter({it[0].pheno_chr== 1.0d})

    pb_quasar_gwas = RUN_QUASAR_PB_GWAS(pb_quasar_gwas_input)

    sc_quasar_gwas_file = Channel
        .of("dataset\tcell_type\tgeno_chr\tpheno_chr\tsig_variant_file\ttime_file\tn_variants")
        .concat(sc_quasar_gwas
            .map { info, sig_variant_file, time_file, n_variants_file ->
                "${info.dataset}\t${info.cell_type}\t${info.chr}\tchr${info.pheno_chr}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(name: 'sc_quasar_gwas_file', newLine: true, sort: false)

    pb_quasar_gwas_file = Channel
        .of("dataset\tcell_type\tmodel\tgeno_chr\tpheno_chr\tsig_variant_file\ttime_file\tn_variants")
        .concat(pb_quasar_gwas
            .map { info, sig_variant_file, time_file, n_variants_file ->
                "${info.dataset}\t${info.cell_type}\t${info.model}\t${info.chr}\tchr${info.pheno_chr}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(name: 'pb_quasar_gwas_file', newLine: true, sort: false)

    csaqtl_quasar_file = Channel
        .of("dataset\tcell_type\tgrouping\tgeno_chr\tsig_variant_file\ttime_file\tn_variants")
        .concat(csaqtl_quasar
            .map { info, sig_variant_file, time_file, n_variants_file ->
                "${info.dataset}\t${info.cell_type}\t${info.grouping}\t${info.chr}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(name: 'csaqtl_quasar_file',newLine: true, sort: false)

    pc_gwas_file = Channel
        .of("dataset\tcell_type\tmodel\tsig_variant_file\ttime_file\tn_variants")
        .concat(pc_gwas
            .map { info, sig_variant_file, time_file, n_variants_file ->
                "${info.dataset}\t${info.cell_type}\t${info.model}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(name: 'pc_gwas_file', newLine: true, sort: false)

    pc_sc_gwas_file = Channel
        .of("dataset\tcell_type\tmodel\tsig_variant_file\ttime_file\tn_variants")
        .concat(pc_sc_gwas
            .map { info, sig_variant_file, time_file, n_variants_file ->
                "${info.dataset}\t${info.cell_type}\t${info.model}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(name: 'pc_sc_gwas_file', newLine: true, sort: false)

    // Genotype permutation analysis.
    rep_bed_files = bed_files
      .combine(channel.of(1..10))
    permute_bed_files = PERMUTE_BED(rep_bed_files)

    perm_sc_quasar_input = sc_quasar_input_full
      .filter( {it[0].data_type != "sct_counts"} )
      .filter { it[0].cov_spec == "bulk_pca" || int_cov_for_cov_spec.containsKey(it[0].cov_spec) }
      .filter( { is_int_sc_cov(it[0].int_cov) } )
      .filter( {it[0].k in ["none", 5]})
      .filter { it[0].cell_type == "T_all" && it[0].int_cov in ["starcat_all", "starcat_CD4_Naive"] }
      // No permutations on the CD4_NC joint cell × individual grid (keep full data only).
      .filter( { !is_joint_both_frac(it[0]) || it[0].is_full_data } )
      .map { info, sc_pheno, covs, anno, bed, cg ->
        tuple([info.dataset, info.chr], info, sc_pheno, covs, anno, bed, cg)
       }
      .combine(permute_bed_files.map { ds, chr, pb -> tuple([ds, chr], pb) }, by: 0)
      .map { key, info, sc_pheno, covs, anno, bed, cg, perm_bed ->
        tuple(info, sc_pheno, covs, anno, perm_bed, cg)
      }
    
    perm_sc_quasar = RUN_QUASAR_SC_PERM(perm_sc_quasar_input)
    perm_sc_quasar_filt = FILTER_VARIANTS_SC(Utils.combineWithPrunedSnps(perm_sc_quasar, pruned_snps))
    perm_sc_quasar_filt = Utils.attachGeneProperties(perm_sc_quasar_filt, gene_properties)

    rep_bed_files_pb = bed_files
      .combine(channel.of(1..5))
    permute_bed_files_pb = PERMUTE_BED_PB(rep_bed_files_pb)

    perm_pb_quasar_input = pb_quasar_input
      .filter( {it[0].cell_type in ["Plasma", "B_IN"]})
      .filter( {it[0].is_full_data} )
      // No permutations on the CD4_NC joint cell × individual grid (keep full data only).
      .filter( { !is_joint_both_frac(it[0]) || it[0].is_full_data } )
      .map { info, pheno, covs, bed, grm ->
        tuple([info.dataset, info.chr], info, pheno, covs, bed, grm)
       }
      .combine(permute_bed_files_pb.map { ds, chr, pb -> tuple([ds, chr], pb) }, by: 0)
      .map { key, info, pheno, covs, bed, grm, perm_bed ->
        tuple(info, pheno, covs, perm_bed, grm)
      }

    perm_pb_quasar = RUN_QUASAR_PB_PERM(perm_pb_quasar_input)
    perm_pb_quasar_filt = FILTER_VARIANTS_PB(Utils.combineWithPrunedSnps(perm_pb_quasar, pruned_snps))
    perm_pb_quasar_filt = Utils.attachGeneProperties(perm_pb_quasar_filt, gene_properties)

    // GWAS permutations.
    gwas_rep_bed_files = bed_files
        .combine(Channel.value(1))
    gwas_permute_bed_files = PERMUTE_BED_GWAS(gwas_rep_bed_files)

    pb_quasar_gwas_perm_input = pb_quasar_input
        .filter({it[0].is_full_data})
        .filter({it[0].cell_type == "B_IN"})
        .filter({it[0].int_cov == "none"})
        .filter({it[0].model == "nb_glm"})
        .combine(Channel.from(1..22))
        .map { info, pheno, covs, bed, grm, pheno_chr ->
          tuple(info + [pheno_chr: pheno_chr as int], pheno, covs, bed)
        }
        .map { info, pheno, covs, bed ->
          tuple([info.dataset, info.chr], info, pheno, covs, bed)
        }
        .combine(gwas_permute_bed_files.map { ds, chr, pb -> tuple([ds, chr], pb) }, by: 0)
        .map({ it -> [it[1], it[2], it[3], it[5]]})
        .filter({it[0].pheno_chr== 1.0d})

    pb_quasar_gwas_perm = RUN_QUASAR_PB_GWAS_PERM(pb_quasar_gwas_perm_input)

    // Disabled for now.
    // sc_quasar_gwas_perm_input = sc_quasar_input
    //     .filter({it[0].model == "lmm_sc"})
    //     .filter({it[0].cov_spec == "bulk_pca"})
    //     .filter({it[0].is_full_data})
    //     .filter({it[0].cell_type == "B_IN"})
    //     .combine(Channel.from(1..22))
    //     .map { info, sc_pheno, covs, anno, bed, pheno_chr ->
    //      tuple(info + [pheno_chr: pheno_chr as int], sc_pheno, covs, anno, bed)
    //     }
    //     .map { info, sc_pheno, covs, anno, bed ->
    //       tuple([info.dataset, info.chr], info, sc_pheno, covs, anno, bed)
    //     }
    //     .combine(gwas_permute_bed_files.map { ds, chr, pb -> tuple([ds, chr], pb) }, by: 0)
    //     .map({ it -> [it[1], it[2], it[3], it[4], it[6]]})
    //     .filter({it[0].pheno_chr== 1.0d})
    // sc_quasar_gwas_perm = RUN_QUASAR_SC_GWAS_PERM(sc_quasar_gwas_perm_input)
    sc_quasar_gwas_perm = Channel.empty()

    sc_quasar_gwas_perm_file = Channel
        .of("dataset\tcell_type\tgeno_chr\tpheno_chr\tsig_variant_file\ttime_file\tn_variants")
        .concat(sc_quasar_gwas_perm
            .map { info, sig_variant_file, time_file, n_variants_file ->
                "${info.dataset}\t${info.cell_type}\t${info.chr}\tchr${info.pheno_chr}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(name: 'sc_quasar_gwas_perm_file', newLine: true, sort: false)

    pb_quasar_gwas_perm_file = Channel
        .of("dataset\tcell_type\tmodel\tgeno_chr\tpheno_chr\tsig_variant_file\ttime_file\tn_variants")
        .concat(pb_quasar_gwas_perm
            .map { info, sig_variant_file, time_file, n_variants_file ->
                "${info.dadtaset}\t${info.cell_type}\t${info.model}\t${info.chr}\tchr${info.pheno_chr}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(name: 'pb_quasar_gwas_perm_file', newLine: true, sort: false)

    // Permutations under an interaction null.
    rep_pb_cov_files = pb_covs
      .filter({it[0].is_full_data})
      .combine(int_cov.filter({it != "none"}))
      .combine(channel.of(1..5))
      .map { info, covs, int_cov, ind ->
        tuple(info.dataset, int_cov, info.cell_type, covs, ind)
      }

    permute_pb_covs = PERMUTE_COV(rep_pb_cov_files)

    perm_int_pb_quasar_input = pb_quasar_input
       .filter( {it[0].is_full_data} )
       .filter( {it[0].cell_type in ["Plasma", "B_IN"]} )
       .filter( {it[0].model in ["nb_glm", "nb_glmm", "lmm", "lm", "p_glmm"]} )
       .map { info, pheno, covs, bed, grm ->
        tuple(info.dataset, info.int_cov, info.cell_type, info, pheno, covs, bed, grm)
       }
       .combine(permute_pb_covs, by: [0, 1, 2])
       .map { dataset, int_cov, cell_type, info, pheno, covs, bed, grm, perm_cov ->
        tuple(info + [interaction_cov: "${info.int_cov}_perm"], pheno, perm_cov, bed, grm)
       }

    perm_int_pb_quasar = RUN_QUASAR_PB_PERM_INT(perm_int_pb_quasar_input)
    perm_int_pb_quasar_filt = FILTER_VARIANTS_PB_INT(Utils.combineWithPrunedSnps(perm_int_pb_quasar, pruned_snps))
    perm_int_pb_quasar_filt = Utils.attachGeneProperties(perm_int_pb_quasar_filt, gene_properties)

    perm_int_pb_quasar_file = Channel
        .of("dataset\tcell_type\tmodel\tchr\tcell_frac\tindiv_frac\tint_cov\tpb_type\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_int_pb_quasar_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.dataset}\t${info.cell_type}\t${info.model}\t${info.chr}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.pb_type}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(name: 'pb_quasar_perm_int_file', newLine: true, sort: false)

    rep_sc_cov_files = sc_covs
      .filter { info, cov_spec, cov ->
          int_cov_for_cov_spec.containsKey(cov_spec) && !is_starcat_int_cov(int_cov_for_cov_spec[cov_spec])
      }
      .map { info, cov_spec, cov -> tuple(info.dataset, int_cov_for_cov_spec[cov_spec], info.cell_type, cov) }
      .combine(channel.of(1..10))
      .map { dataset, int_cov, cell_type, covs, ind -> tuple(dataset, int_cov, cell_type, covs, ind) }

    permute_sc_covs = PERMUTE_SC_COV(rep_sc_cov_files)

    perm_int_sc_quasar_input = sc_quasar_input_full
      .filter { it[0].data_type != "sct_counts" }
      .filter { it[0].model == "p_glmm_sc" }
      .filter { it[0].k == "none" }
      .filter { it[0].int_cov != "none" && !is_starcat_int_cov(it[0].int_cov) }
      .map { info, sc_pheno, covs, anno, bed, cg ->
        tuple(info.dataset, info.int_cov, info.cell_type, info, sc_pheno, covs, anno, bed, cg)
      }
      .combine(permute_sc_covs, by: [0, 1, 2])
      .map { dataset, int_cov, cell_type, info, sc_pheno, covs, anno, bed, cg, perm_cov ->
       tuple(info + [interaction_cov: "${info.int_cov}_perm"], sc_pheno, perm_cov, anno, bed, cg)
      }

    perm_int_sc_quasar = RUN_QUASAR_SC_PERM_INT(perm_int_sc_quasar_input)
    perm_int_sc_quasar_filt = FILTER_VARIANTS_SC_INT(Utils.combineWithPrunedSnps(perm_int_sc_quasar, pruned_snps))
    perm_int_sc_quasar_filt = Utils.attachGeneProperties(perm_int_sc_quasar_filt, gene_properties)

    sc_quasar_perm_int_file = Channel
        .of("dataset\tmodel\tcell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_int_sc_quasar_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.dataset}\t${info.model}\t${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(name: 'sc_quasar_perm_int_file', newLine: true, sort: false)

    rep_sc_cov_files_within = sc_covs
      .filter { info, cov_spec, cov ->
          int_cov_for_cov_spec.containsKey(cov_spec) && !is_starcat_int_cov(int_cov_for_cov_spec[cov_spec])
      }
      .map { info, cov_spec, cov -> tuple(info.dataset, int_cov_for_cov_spec[cov_spec], info.cell_type, cov) }
      .combine(channel.of(1..10))
      .map { dataset, int_cov, cell_type, covs, ind -> tuple(dataset, int_cov, cell_type, covs, ind) }

    permute_sc_covs_within = PERMUTE_SC_COV_WITHIN(rep_sc_cov_files_within)

    perm_int_sc_quasar_within_input = sc_quasar_input_full
      .filter({ it[0].data_type != "sct_counts"} )
      .filter({ it[0].model == "p_glmm_sc"} )
      .filter({ it[0].k == "none"} )
      .filter({ it[0].int_cov != "none" && !is_starcat_int_cov(it[0].int_cov)})
      .map { info, sc_pheno, covs, anno, bed, cg ->
        tuple(info.dataset, info.int_cov, info.cell_type, info, sc_pheno, covs, anno, bed, cg)
      }
      .combine(permute_sc_covs_within, by: [0, 1, 2])
      .map { dataset, int_cov, cell_type, info, sc_pheno, covs, anno, bed, cg, perm_cov ->
        tuple(info + [interaction_cov: "${info.int_cov}_perm"], sc_pheno, perm_cov, anno, bed, cg)
      }

    perm_int_sc_quasar_within = RUN_QUASAR_SC_PERM_WITHIN_INT(perm_int_sc_quasar_within_input)
    perm_int_sc_quasar_within_filt = FILTER_VARIANTS_SC_INT_WITHIN(Utils.combineWithPrunedSnps(perm_int_sc_quasar_within, pruned_snps))
    perm_int_sc_quasar_within_filt = Utils.attachGeneProperties(perm_int_sc_quasar_within_filt, gene_properties)

    sc_quasar_perm_int_within_file = Channel
        .of("dataset\tmodel\tcell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_int_sc_quasar_within_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.dataset}\t${info.model}\t${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(name: 'sc_quasar_perm_int_within_file', newLine: true, sort: false)

    // Shuffle the interaction covariate within log-library-size bins so that
    // depth-associated pseudotime structure is preserved.
    rep_sc_cov_files_logcount = sc_covs
      .filter { info, cov_spec, cov ->
          int_cov_for_cov_spec.containsKey(cov_spec) && !is_starcat_int_cov(int_cov_for_cov_spec[cov_spec])
      }
      .map { info, cov_spec, cov ->
        tuple(
          [info.dataset, info.cell_type, info.cell_frac, info.indiv_frac, info.n_cells_target, info.count_frac],
          info.dataset,
          int_cov_for_cov_spec[cov_spec],
          info.cell_type,
          cov
        )
      }
      .combine(
        sc_counts.map { type, info, counts ->
          tuple(
            [info.dataset, info.cell_type, info.cell_frac, info.indiv_frac, info.n_cells_target, info.count_frac],
            counts
          )
        },
        by: 0
      )
      .combine(channel.of(1..10))
      .map { key, dataset, int_cov, cell_type, cov, counts, ind ->
        tuple(dataset, int_cov, cell_type, cov, counts, ind)
      }

    permute_sc_covs_logcount = PERMUTE_SC_COV_LOGCOUNT_BINS(rep_sc_cov_files_logcount)

    perm_int_sc_quasar_logcount_input = sc_quasar_input_full
      .filter({ it[0].data_type in ["counts", "sct_counts"]} )
      .filter({ it[0].model == "p_glmm_sc"} )
      .filter({ it[0].k == "none"} )
      .filter({ it[0].int_cov != "none" && !is_starcat_int_cov(it[0].int_cov)})
      // TEMP: skip interaction perms while only running the unpermuted quasar interaction.
      .filter { false }
      .map { info, sc_pheno, covs, anno, bed, cg ->
        tuple(info.dataset, info.int_cov, info.cell_type, info, sc_pheno, covs, anno, bed, cg)
      }
      .combine(permute_sc_covs_logcount, by: [0, 1, 2])
      .map { dataset, int_cov, cell_type, info, sc_pheno, covs, anno, bed, cg, perm_cov ->
        tuple(info + [interaction_cov: "${info.int_cov}_perm"], sc_pheno, perm_cov, anno, bed, cg)
      }

    perm_int_sc_quasar_logcount = RUN_QUASAR_SC_PERM_LOGCOUNT_INT(perm_int_sc_quasar_logcount_input)
    perm_int_sc_quasar_logcount_filt = FILTER_VARIANTS_SC_INT_LOGCOUNT(Utils.combineWithPrunedSnps(perm_int_sc_quasar_logcount, pruned_snps))
    perm_int_sc_quasar_logcount_filt = Utils.attachGeneProperties(perm_int_sc_quasar_logcount_filt, gene_properties)

    sc_quasar_perm_int_logcount_file = Channel
        .of("dataset\tmodel\tdata_type\tcell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_int_sc_quasar_logcount_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.dataset}\t${info.model}\t${info.data_type}\t${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(name: 'sc_quasar_perm_int_logcount_file', newLine: true, sort: false)

    perm_sc_quasar_file = Channel
        .of("dataset\tmodel\tcell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tcount_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_sc_quasar_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.dataset}\t${info.model}\t${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.count_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(name: 'sc_quasar_perm_file', newLine: true, sort: false)

    perm_pb_quasar_file = Channel
        .of("dataset\tcell_type\tmodel\tchr\tcell_frac\tindiv_frac\tcount_frac\tint_cov\tpb_type\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_pb_quasar_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.dataset}\t${info.cell_type}\t${info.model}\t${info.chr}\t${info.cell_frac}\t${info.indiv_frac}\t${info.count_frac}\t${info.int_cov}\t${info.pb_type}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(name: 'pb_quasar_perm_file', newLine: true, sort: false)

    // Analysis.
    //PLOT_CLUMPED(clumped_pb_quasar_file, clumped_sc_quasar_file)
    //PLOT_POWER(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    //PLOT_POWER_N_CELLS_FILTER(pb_quasar_file, sc_quasar_file)
    //PLOT_POWER_OFFSET(sc_offset_quasar_file)
    //unique_egenes_tsv = PLOT_UNIQUE_SC_EGENES(pb_quasar_file, sc_quasar_file)
    //PLOT_EGENE_SIG_MODEL(pb_quasar_file, sc_quasar_file)
    //PLOT_POWER_CELL_FRAC(pb_quasar_file, sc_quasar_file)
    // TEMP: disabled with joint downsampling
    //PLOT_POWER_BOTH_FRAC(pb_quasar_file, sc_quasar_file)
    //PLOT_POWER_COUNT_FRAC(pb_quasar_file, sc_quasar_file)
    //PLOT_POWER_INDIV_FRAC(pb_quasar_file, sc_quasar_file)
    //PLOT_POWER_GENE_PROP(pb_quasar_file, sc_quasar_file)
    //PLOT_METHOD_SCATTER(pb_quasar_file, sc_quasar_file)

    //unique_genotype_dosages = EXTRACT_GENOTYPES(
    //    all_bed.combine(unique_egenes_tsv)
    //)
    //unique_sc_logcounts = EXTRACT_SC_LOGCOUNTS(
    //   sc_preprocess_input
    //        .filter { it[0].is_full_data && it[0].cell_type in ["Plasma", "B_IN", "CD4_NC"] }
    //       .combine(unique_egenes_tsv)
    //)
    //unique_sc_egene_figures_input = unique_sc_logcounts
    //    .map { info, logcounts -> tuple(info.dataset, info, logcounts) }
    //    .combine(unique_genotype_dosages, by: 0)
    //    .combine(unique_egenes_tsv)
    //    .map { dataset, info, logcounts, geno, leads ->
    //        tuple(info, leads, logcounts, geno)
    //    }
    //PLOT_UNIQUE_SC_EGENE_FIGURES(unique_sc_egene_figures_input)

    //unique_pb_egenes_tsv = PLOT_UNIQUE_PB_EGENES(pb_quasar_file, sc_quasar_file)
    //unique_pb_genotype_dosages = EXTRACT_GENOTYPES_PB_UNIQUE(
    //    all_bed.combine(unique_pb_egenes_tsv)
    //)
    //unique_pb_sc_logcounts = EXTRACT_SC_LOGCOUNTS_PB_UNIQUE(
    //    sc_preprocess_input
    //        .filter { it[0].is_full_data && it[0].cell_type in ["Plasma", "B_IN", "CD4_NC"] }
    //        .combine(unique_pb_egenes_tsv)
    //)
    //unique_pb_egene_figures_input = unique_pb_sc_logcounts
    //    .map { info, logcounts -> tuple(info.dataset, info, logcounts) }
    //    .combine(unique_pb_genotype_dosages, by: 0)
    //    .combine(unique_pb_egenes_tsv)
    //    .map { dataset, info, logcounts, geno, leads ->
    //        tuple(info, leads, logcounts, geno)
    //    }
    //PLOT_UNIQUE_PB_EGENE_FIGURES(unique_pb_egene_figures_input)
    //PLOT_POWER_COVS(sc_quasar_file)
    //PLOT_CONVERGENCE(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    //PLOT_PERM(perm_sc_quasar_file)
    //PLOT_PERM_CELL_FRAC(perm_sc_quasar_file, perm_pb_quasar_file)
    //PLOT_PERM_COUNT_FRAC(perm_sc_quasar_file, perm_pb_quasar_file)
    //PLOT_PERM_PB(perm_pb_quasar_file)
    PLOT_PERM_GLOBAL(perm_pb_quasar_file)
    //PLOT_PERM_GWAS(sc_quasar_gwas_perm_file, pb_quasar_gwas_perm_file)
    PLOT_PERM_INT(perm_int_pb_quasar_file)
    PLOT_INT_OUTPUT(pb_quasar_file)
    //PLOT_SC_INT_OUTPUT(sc_quasar_file, castie_file)
    PLOT_QUASAR_SC_INT_OUTPUT(sc_quasar_file)
    //PLOT_GROUPED_SC_OUTPUT(sc_quasar_file)
    //PLOT_GROUPED_INT_TIME(sc_quasar_file)
    //PLOT_GROUPED_VS_INT(sc_quasar_file, castie_file)
    //PLOT_SC_PGLMM_VS_LMM(sc_quasar_file)
    //PLOT_PERM_SC_GROUPED(perm_sc_quasar_file)
    PLOT_PERM_SC_INT_GLOBAL(perm_sc_quasar_file)
    //PLOT_PERM_SC_INT(sc_quasar_perm_int_file)
    //PLOT_SC_INT_FIGURES(sc_int_figures_input)
    //PLOT_METACELL_OUTPUT(sc_quasar_file, seacells_file)
    //PLOT_CSAQTL_OUTPUT(csaqtl_quasar_file)
    //PLOT_PC_GWAS_OUTPUT(pc_gwas_file, pc_sc_gwas_file)so
    //PLOT_PERM_SC_INT_WITHIN(sc_quasar_perm_int_within_file)
    //PLOT_PERM_SC_INT_LOGCOUNT_BINS(sc_quasar_perm_int_logcount_file)
    //PLOT_TIME(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    //PLOT_CELLS_PER_INDIV(
    //    cluster_sizes.map { dataset, sizes, cells_per_indiv -> cells_per_indiv }
    //)
    //PLOT_PC_GWAS_COMPARISON(pc_gwas_file, pc_sc_gwas_file)
    //PLOT_CONCORDANCE(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    //PLOT_GENE_PROPERTIES(gene_properties_file)
    //PLOT_GWAS_OUTPUT(sc_quasar_gwas_file, pb_quasar_gwas_file)
    //PLOT_PVALUE_SCATTER(pb_quasar_file, sc_quasar_file)
    //PLOT_ZSCORE_SCATTER(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    //PLOT_MAIN_VS_INT(sc_quasar_file)
    //PLOT_INT_TIME(sc_quasar_file, castie_file)
    PLOT_QUASAR_INT_TIME(sc_quasar_file)
    //RUN_COLOC(sc_quasar_file)

    // Create example data for quasar.
    //example_data_input = sc_quasar_input
    //    .filter({it[0].cell_type == "B_IN"})
    //    .filter({it[0].indiv_frac == 1.0d && it[0].cell_frac == 1.0d})
    //    .filter({it[0].chr == "chr22"})
    //CREATE_EXAMPLE_DATA(example_data_input)
}
