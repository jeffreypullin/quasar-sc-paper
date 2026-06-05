#!/usr/bin/env nextflow

params.onek1k_raw_single_cell_data = file("data/onek1k/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz")
params.onek1k_supp_tables = file("data/onek1k/science.abf3041_tables_s6_to_s19/science.abf3041_tables_s6_to_s19.xlsx")
params.annot_bed_url = "https://ftp.ensembl.org/pub/grch37/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
params.onek1k_cell_types = [
    'CD4 NC', 'CD4 ET', 'CD4 SOX4', 'CD8 ET', 'CD8 NC', 'CD8 S100B', 'NK', 
    'NK R', 'Plasma', 'B Mem', 'B IN', 'Mono C', 'Mono NC', 'DC'
]

include { 
    EXTRACT_INDIV_IDS ; FILTER_VCF ; CONVERT_VCF_TO_BED ; PRUNE_SNPS ; CONCAT_BED_FILES ;
    EXTRACT_GENOTYPES ; COMPUTE_SC_COUNTS ; EXTRACT_SC_LOGCOUNTS ; COMPUTE_PB_COUNTS ; COMPUTE_CLUSTER_SIZES ; COMPUTE_PB_EXPR_COVS ;
    COMPUTE_SC_EXPR_COVS ; COMPUTE_GENOTYPE_PCS ; CREATE_ANNOT_BED ; COMPUTE_PB_LOGCOUNTS ;
    CREATE_GRM ; PREPARE_SLINGSHOT_ADATA ; PREPARE_SEACELLS_ADATA ; RUN_SLINGSHOT ; JOIN_SC_INT_COVS ; JOIN_PB_INT_COVS ;
    COMPUTE_STARCAT_COVS ; COMPUTE_CSAQTL_COUNTS
} from './modules/preprocess'
include { COLLATE_SAIGEQTL_INPUT ; SAIGE_SUBSET_BED ; EXTRACT_GENES ; 
          RUN_SAIGEQTL ; CONCAT_SAIGEQTL_TSVS ; COMPUTE_SAIGEQTL_POWER } from './modules/saige-qtl'
include { COMBINE_PB_COVS ; COMBINE_SC_COVS ; FILTER_COVS ; ANNOTATE_PHENO ; 
          RUN_QUASAR_PB ; RUN_QUASAR_SC ; FILTER_VARIANTS ; RUN_QUASAR_SC_GWAS ;
          RUN_QUASAR_PB_GWAS ; COMPUTE_CELL_GROUPS ; COMPUTE_SEACELLS ;
          RUN_CSAQTL_QUASAR } from './modules/quasar'
include { COMPUTE_CELL_GROUPS as COMPUTE_PSEUDOTIME_GROUPS } from './modules/quasar'
include { FILTER_VARIANTS as FILTER_VARIANTS_SC } from './modules/quasar'
include { FILTER_VARIANTS as FILTER_VARIANTS_PB } from './modules/quasar'
include { FILTER_VARIANTS as FILTER_VARIANTS_PB_INT } from './modules/quasar'
include { FILTER_VARIANTS as FILTER_VARIANTS_SC_INT } from './modules/quasar'
include { FILTER_VARIANTS as FILTER_VARIANTS_SC_INT_WITHIN } from './modules/quasar'
include { COMPUTE_QUASAR_POWER as COMPUTE_QUASAR_POWER_SC } from './modules/quasar'
include { COMPUTE_QUASAR_POWER as COMPUTE_QUASAR_POWER_PB } from './modules/quasar'
include { RUN_QUASAR_PB as RUN_QUASAR_PB_PERM } from './modules/quasar'
include { RUN_QUASAR_SC as RUN_QUASAR_SC_PERM } from './modules/quasar'
include { RUN_QUASAR_PB as RUN_QUASAR_PB_PERM_INT } from './modules/quasar'
include { RUN_QUASAR_SC as RUN_QUASAR_SC_PERM_INT } from './modules/quasar'
include { RUN_QUASAR_SC as RUN_QUASAR_SC_PERM_WITHIN_INT } from './modules/quasar'
include { RUN_QUASAR_PB_GWAS as RUN_QUASAR_PB_GWAS_PERM } from './modules/quasar'
include { RUN_QUASAR_SC_GWAS as RUN_QUASAR_SC_GWAS_PERM } from './modules/quasar'
include { PERMUTE_COV as PERMUTE_SC_COV } from './modules/analysis'
include { 
    PERMUTE_BED ; PLOT_POWER ; PLOT_CONVERGENCE ; 
    PLOT_PERM ; PLOT_PERM_PB ; PLOT_TIME ; PLOT_CONCORDANCE ;
    COMPUTE_GENE_PROPERTIES ; PLOT_GENE_PROPERTIES ; PLOT_PERM_INT ;
    PLOT_INT_OUTPUT ; PLOT_GWAS_OUTPUT ; CREATE_EXAMPLE_DATA ;
    PLOT_PERM_GLOBAL ; PERMUTE_COV ; PLOT_PERM_GWAS ;
    PLOT_PVALUE_SCATTER ; PLOT_SC_INT_OUTPUT ; PLOT_PERM_SC_INT ; PLOT_PERM_SC_INT_WITHIN ;
    //CLUMP_VARIANTS ;
    PLOT_PERM_SC_GROUPED ; PLOT_GROUPED_SC_OUTPUT ; PLOT_PERM_SC_INT_GLOBAL ;
    //PLOT_CLUMPED ;
    PERMUTE_SC_COV_WITHIN ; PLOT_SC_INT_FIGURES ; PLOT_METACELL_OUTPUT ; PLOT_CSAQTL_OUTPUT
    } from './modules/analysis'
//include { CLUMP_VARIANTS as CLUMP_VARIANTS_PB } from "./modules/analysis"
//include { CLUMP_VARIANTS as CLUMP_VARIANTS_SC } from "./modules/analysis"
include { PERMUTE_BED as PERMUTE_BED_GWAS} from './modules/analysis'
include { COMPUTE_CONVERGENCE as COMPUTE_CONVERGENCE_PB } from './modules/analysis'
include { COMPUTE_CONVERGENCE as COMPUTE_CONVERGENCE_SC } from './modules/analysis'

workflow {

    vcf_files = channel.fromFilePairs(
        "data/onek1k/chr*.dose.filtered.R2_0.8.vcf.gz",
        size: 1, 
        flat: true
    )

    ids = EXTRACT_INDIV_IDS(params.onek1k_raw_single_cell_data)
    cluster_sizes = COMPUTE_CLUSTER_SIZES(params.onek1k_raw_single_cell_data)
    filt_vcf_files = FILTER_VCF(vcf_files, ids)
    bed_files = CONVERT_VCF_TO_BED(filt_vcf_files)
    pruned_snps = PRUNE_SNPS(bed_files)
    all_bed = CONCAT_BED_FILES(bed_files.map({ it[1] }).collect())
    geno_pcs = COMPUTE_GENOTYPE_PCS(all_bed)
    grm = CREATE_GRM(all_bed)
    anno_bed = CREATE_ANNOT_BED(params.annot_bed_url)

    cell_infos = channel.fromList(params.onek1k_cell_types)
        .flatMap { ct ->
            def info = [cell_type: ct.toString(), cell_frac: 1.0d, indiv_frac: 1.0d]
            def combos = []
            if (info.cell_type in ["Plasma", "B IN", "CD4 NC"]) {
                def cell_fracs = [0.005d, 0.01d, 0.05d, 0.1d, 0.25d, 0.5d, 1.0d]
                combos += cell_fracs.collect { f -> info + [cell_frac: f, indiv_frac: 1.0d] }
                if (info.cell_type == "CD4 NC") {
                    def indiv_fracs = [0.05d, 0.1d, 0.25d, 0.5d, 1.0d]
                    combos += indiv_fracs.collect { f -> info + [cell_frac: 1.0d, indiv_frac: f] }
                    return combos.unique { x -> "${x.cell_frac}|${x.indiv_frac}" }
                }
                return combos
            }
            return [info]
        }

    cell_infos = cell_infos.mix(
        channel.fromList([
            [cell_type: 'B_all', cell_frac: 1.0d, indiv_frac: 1.0d],
            [cell_type: 'CD4_T_all', cell_frac: 1.0d, indiv_frac: 1.0d],
            [cell_type: 'CD8_T_all', cell_frac: 1.0d, indiv_frac: 1.0d],
        ])
    )

    sc_input = cell_infos
        .combine(Channel.of(params.onek1k_raw_single_cell_data))

    gene_properties = COMPUTE_GENE_PROPERTIES(sc_input, anno_bed)
    gene_properties_plot = gene_properties.filter { info, properties_file ->
        info.cell_frac == 1.0d && info.indiv_frac == 1.0d
    }
    gene_properties_file = Channel
        .of("cell_type\tproperties_file")
        .concat(gene_properties_plot
            .map { info, properties_file->
                "${info.cell_type}\t${properties_file}"
            }
        )
        .collectFile(
            name: 'gene_properties_file',
            newLine: true,
            sort: false,
            storeDir: "${projectDir}/data/intermediate"
        )

    sc_input = sc_input.combine(gene_properties, by: 0)
    pb_counts = COMPUTE_PB_COUNTS(sc_input)
    pb_logcounts = COMPUTE_PB_LOGCOUNTS(sc_input)
    pb_input = pb_counts.mix(pb_logcounts)
    pheno = ANNOTATE_PHENO(pb_input, anno_bed)
    sc_counts = COMPUTE_SC_COUNTS(sc_input)

    pb_expr_covs = COMPUTE_PB_EXPR_COVS(sc_input)
    sc_expr_covs = COMPUTE_SC_EXPR_COVS(sc_input)
    starcat_covs = COMPUTE_STARCAT_COVS(
        sc_input.filter {
            it[0].cell_type in ["CD4_T_all", "CD8_T_all"] &&
            it[0].cell_frac == 1.0d &&
            it[0].indiv_frac == 1.0d
        }
    )

    sc_input = sc_counts
        .combine(pb_expr_covs, by: 0)
        .combine(geno_pcs)

    // Run SAIGE-QTL.
    saige_subset_bed = SAIGE_SUBSET_BED(all_bed)
    saigeqtl_input = sc_input
        .combine(bed_files)
        .combine(anno_bed)
        .map({ x ->
            def dict = x[0] + [ chr: x[4] ]
            tuple(dict, x[1], x[2], x[3], x[5], x[6])
        })
        .filter( { it[0].cell_type in ["Plasma"]} )
        | COLLATE_SAIGEQTL_INPUT

    gene_lists = EXTRACT_GENES(saigeqtl_input.map({[it[0], it[1]]}), 50)
        .flatMap { info, gene_lists ->
            def lists = (gene_lists instanceof List) ? gene_lists : [gene_lists]
            lists.collect { glist -> tuple(info, glist) }
        }

    saigeqtl_input = saigeqtl_input 
        .combine(gene_lists, by: 0)
        .filter { it[0].cell_frac == 1.0d && it[0].indiv_frac == 1.0d }
        .combine(saige_subset_bed)
        .map { info, input, bed, annot, gene_file, subset_bed ->
            def new_info = info + [cell_type: info.cell_type.replaceAll(/\s+/, '-')]
            tuple(new_info, input, bed, annot, gene_file, subset_bed)
        }

    saigeqtl_output = RUN_SAIGEQTL(saigeqtl_input)
    saigeqtl_output = COMPUTE_SAIGEQTL_POWER(saigeqtl_output)
    saigeqtl_file = CONCAT_SAIGEQTL_TSVS(saigeqtl_output.map{it[1]}.collect())

    // Run pseudobulk quasar.
    pb_covs = COMBINE_PB_COVS(pb_expr_covs, geno_pcs)
    model = Channel.of("nb_glm", "lm", "lmm")
    int_cov = Channel.of("none", "age", "sex")
    pb_type = Channel.of("counts", "logcounts")
    pb_quasar_input = pheno
        .combine(pb_covs, by: 0)
        .combine(bed_files)
        .combine(model)
        .combine(int_cov)
        .combine(grm) 
        .map({ x ->
            def dict = x[0] + [ chr: x[4], model: x[6], int_cov: x[7], pb_type: x[1]]
            tuple(dict, x[2], x[3], x[5], x[8])
        })
        .filter({(!(it[0].model in ["lm", "lmm"]) && it[0].pb_type == "counts") || 
                 ((it[0].model in ["lm", "lmm"]) && it[0].pb_type == "logcounts") })
        .filter({it[0].indiv_frac == 1.0d && it[0].cell_frac == 1.0d })
        .filter({it[0].model in ["lm", "nb_glm", "lmm"]})
        //.filter({it[0].int_cov != "none" && it[0].model == "nb_glm"})
        //.filter({it[0].cell_type in ["Plasma", "B IN", "CD4 NC"]})
        //.filter({it[0].model in ["nb_glm", "lm"]})
        .filter({it[0].cell_type in ["Plasma"]})
        .filter({it[0].model in ["lm"]})

    // Slingshot pseudotime for B_all (used by both sc and pb pseudotime branches).
    slingshot_input = cell_infos
        .filter { it.cell_type == "B_all" }
        .combine(Channel.of(params.onek1k_raw_single_cell_data))
    slingshot_adata = PREPARE_SLINGSHOT_ADATA(slingshot_input)
    slingshot_out = RUN_SLINGSHOT(slingshot_adata)

    // SEACells metacells for all full-fraction clusters (independent of slingshot).
    seacells_input = cell_infos
        .filter { it.cell_frac == 1.0d && it.indiv_frac == 1.0d }
        .filter { !(it.cell_type in ["Plasma", "CD4 SOX4"]) }
        .combine(Channel.of(params.onek1k_raw_single_cell_data))
    seacells_adata = PREPARE_SEACELLS_ADATA(seacells_input)
    seacells_out = COMPUTE_SEACELLS(seacells_adata)
    seacells_groups = seacells_out.groups

    // Cell-state abundance QTL phenotypes are computed later (after join_sc_pseudo_out
    // is available) so that B_all can use pseudotime-based cell groups.

    // Pseudobulk B_all pseudotime: aggregate per-cell pseudotime to per-sample mean
    // and append to the pb cov file.
    pb_covs_b_all = pb_covs.filter {
        it[0].cell_type == "B_all" &&
        it[0].cell_frac == 1.0d &&
        it[0].indiv_frac == 1.0d
    }
    sc_expr_covs_b_all = sc_expr_covs.filter {
        it[0].cell_type == "B_all" &&
        it[0].cell_frac == 1.0d &&
        it[0].indiv_frac == 1.0d
    }
    join_pb_pseudo_input = pb_covs_b_all
        .combine(sc_expr_covs_b_all, by: 0)
        .combine(slingshot_out)
        .map { info, pb_cov, sc_expr_cov, pt, p1, p2 ->
            tuple(info, pb_cov, sc_expr_cov, pt)
        }
    join_pb_pseudo_out = JOIN_PB_INT_COVS(join_pb_pseudo_input)

    pb_quasar_input_b_all_pseudo = pheno
        .filter {
            it[0].cell_type == "B_all" &&
            it[0].cell_frac == 1.0d &&
            it[0].indiv_frac == 1.0d &&
            it[1] == "logcounts"
        }
        .combine(join_pb_pseudo_out, by: 0)
        .combine(bed_files)
        .combine(grm)
        .map { info, pb_type_val, pheno_bed, cov_pt, chr, bed_file, grm_path ->
            def dict = info + [chr: chr, model: "lm", int_cov: "pseudotime", pb_type: "logcounts"]
            tuple(dict, pheno_bed, cov_pt, bed_file, grm_path)
        }

    //pb_quasar_input = pb_quasar_input.mix(pb_quasar_input_b_all_pseudo)
    //pb_quasar_input = pb_quasar_input_b_all_pseudo
    pb_quasar_input = pb_quasar_input

    pb_quasar = RUN_QUASAR_PB(pb_quasar_input)
    pb_quasar_keyed = pb_quasar.map { info, region, variant, time ->
      tuple([cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac], info, region, variant, time)
    }
    pb_quasar = pb_quasar_keyed
        .combine(gene_properties, by: 0)
        .map { base_info, info, region, variant, time, gp_file ->
            tuple(info, region, variant, time, gp_file)
        }

    pb_quasar = COMPUTE_QUASAR_POWER_PB(pb_quasar)
    pb_quasar = COMPUTE_CONVERGENCE_PB(pb_quasar)

    // Run single-cell quasar. 
    both_expr_covs = pb_expr_covs.combine(sc_expr_covs, by: 0)
    sc_covs = COMBINE_SC_COVS(both_expr_covs, geno_pcs)
    cov_spec = Channel.of("bulk_pca", "sc_pca", "bulk_pca+pct_mito", "bulk_pca+cell_cycle")
    int_cov_sc = Channel.of("none")
    sc_covs = FILTER_COVS(sc_covs.combine(cov_spec))

    // Slingshot pseudotime for B_all sc branch (slingshot_out produced above).
    b_all_filtered_covs_for_join = sc_covs.filter {
        it[0].cell_type == "B_all" &&
        it[0].cell_frac == 1.0d &&
        it[0].indiv_frac == 1.0d &&
        it[1] == "bulk_pca"
    }
    join_slingshot_input = b_all_filtered_covs_for_join
    .combine(slingshot_out)
    .map { info, cov_spec, cov_data, int_cov_data, plot1, plot2 ->
        tuple(info, cov_data, int_cov_data)
    }

    // starCAT/TCAT scores for split T cell grouped interaction analyses.
    starcat_int_cov_by_cell_type = [
        "CD4_T_all": "starcat_CD4_Naive",
        "CD8_T_all": "starcat_Cytotoxic",
    ]
    starcat_cell_types = starcat_int_cov_by_cell_type.keySet() as List
    starcat_int_covs = starcat_int_cov_by_cell_type.values() as List
    starcat_filtered_covs_for_join = sc_covs.filter {
        starcat_cell_types.contains(it[0].cell_type) &&
        it[0].cell_frac == 1.0d &&
        it[0].indiv_frac == 1.0d &&
        it[1] == "bulk_pca"
    }
    join_starcat_input = starcat_filtered_covs_for_join
        .combine(starcat_covs, by: 0)
        .map { info, cov_spec, cov_data, starcat_cov_data, starcat_plot ->
            tuple(info, cov_data, starcat_cov_data)
        }

    joined_sc_int_covs = JOIN_SC_INT_COVS(join_slingshot_input.mix(join_starcat_input))
    join_sc_pseudo_out = joined_sc_int_covs.filter { it[0].cell_type == "B_all" }
    join_sc_starcat_out = joined_sc_int_covs.filter { starcat_cell_types.contains(it[0].cell_type) }

    // Cell-state abundance QTL: per-individual cell-group counts -> nb_glm GWAS.
    // For B_all, bin cells by pseudotime rather than SEACells metacells.
    b_all_pseudotime_groups_input = join_sc_pseudo_out
        .filter { it[0].cell_type == "B_all" }
        .map { info, covs ->
            tuple(
                [cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac,
                 int_cov: "pseudotime", k: 5],
                covs
            )
        }
    b_all_pseudotime_groups = COMPUTE_PSEUDOTIME_GROUPS(b_all_pseudotime_groups_input)
        .map { info, cg ->
            tuple([cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac], cg)
        }

    csaqtl_cell_groups = seacells_groups
        .filter { it[0].cell_type != "B_all" }
        .mix(b_all_pseudotime_groups)

    csaqtl_pheno = COMPUTE_CSAQTL_COUNTS(
        csaqtl_cell_groups.combine(seacells_adata, by: 0)
    )

    csaqtl_quasar_input = csaqtl_pheno
        .combine(
            pb_covs.filter {
                it[0].cell_frac == 1.0d &&
                it[0].indiv_frac == 1.0d
            },
            by: 0
        )
        .combine(all_bed)
        .map { info, pheno, covs, bed ->
            tuple(info, pheno, covs, bed)
        }

    csaqtl_quasar = RUN_CSAQTL_QUASAR(csaqtl_quasar_input)

    sc_quasar_input_core = sc_counts
        .combine(sc_covs, by: 0)
        .combine(anno_bed)
        .combine(bed_files)
        .combine(int_cov_sc)
        .map({ x ->
            def dict = x[0] + [ chr: x[5], cov_spec: x[2], int_cov: x[7]]
            tuple(dict, x[1], x[3], x[4], x[6])
        })
        .filter { !(it[0].cov_spec in ["sc_pca", "bulk_pca+pct_mito", "bulk_pca+cell_cycle"])
                || (it[0].cell_type in ["Plasma", "B IN", "CD4 NC"]) }
        .filter({ it[0].cell_frac == 1.0d && it[0].indiv_frac == 1.0d && it[0].cov_spec == "bulk_pca" })

    sc_quasar_input_b_all_pseudo = sc_counts
        .filter {
            it[0].cell_type == "B_all" &&
            it[0].cell_frac == 1.0d &&
            it[0].indiv_frac == 1.0d
        }
        .combine(join_sc_pseudo_out, by: 0)
        .combine(anno_bed)
        .combine(bed_files)
        .map { info, pheno, cov_pt, anno, chr, plink ->
            def dict = info + [chr: chr, cov_spec: "bulk_pca", int_cov: "pseudotime"]
            tuple(dict, pheno, cov_pt, anno, plink)
        }

    sc_quasar_input_starcat = sc_counts
        .filter {
            starcat_cell_types.contains(it[0].cell_type) &&
            it[0].cell_frac == 1.0d &&
            it[0].indiv_frac == 1.0d
        }
        .combine(join_sc_starcat_out, by: 0)
        .combine(anno_bed)
        .combine(bed_files)
        .map { info, pheno, cov_starcat, anno, chr, plink ->
            def dict = info + [
                chr: chr,
                cov_spec: "bulk_pca",
                int_cov: starcat_int_cov_by_cell_type[info.cell_type],
            ]
            tuple(dict, pheno, cov_starcat, anno, plink)
        }

    //sc_quasar_input = sc_quasar_input_core.mix(sc_quasar_input_b_all_pseudo)
    //sc_quasar_input = sc_quasar_input_b_all_pseudo.mix(sc_quasar_input_starcat)
    //sc_quasar_input = sc_quasar_input_core
    //sc_quasar_input = sc_quasar_input_starcat
    //sc_quasar_input = sc_quasar_input_b_all_pseudo
    sc_quasar_input = sc_quasar_input_core.filter { it[0].cell_type == "CD4_T_all" }

    // Cell groups for the single-cell quasar. Three modes, selected by sc_cell_group_k:
    //   - "none"    : no aggregation; quasar runs per-cell (optionally with -i interaction).
    //   - "seacells": precomputed SEACells metacells (independent of any covariate).
    //   - integer k : equal-width bins along an interaction covariate (e.g. pseudotime),
    //                 only valid when an interaction covariate is set.
    sc_cell_group_k = Channel.of("seacells")

    sc_quasar_input_with_k = sc_quasar_input
        .combine(sc_cell_group_k)
        .map { info, sc_pheno, covs, anno, plink, k ->
            tuple(info + [k: k], sc_pheno, covs, anno, plink)
        }

    // No-groups branch: pass covs as a placeholder in the cell_groups slot (unused).
    sc_no_groups_full = sc_quasar_input_with_k
        .filter { it[0].k == "none" }
        .map { info, sc_pheno, covs, anno, plink ->
            tuple(info, sc_pheno, covs, anno, plink, covs)
        }

    // SEACells branch: join each cell type to its precomputed metacell groups.
    sc_seacells_full = sc_quasar_input_with_k
        .filter { it[0].k == "seacells" }
        .map { info, sc_pheno, covs, anno, plink ->
            def base = [cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac]
            tuple(base, info, sc_pheno, covs, anno, plink)
        }
        .combine(seacells_groups, by: 0)
        .map { base, info, sc_pheno, covs, anno, plink, cg ->
            tuple(info, sc_pheno, covs, anno, plink, cg)
        }

    // Interaction-covariate binning branch: only when an interaction covariate exists.
    sc_binned_in = sc_quasar_input_with_k
        .filter { it[0].k != "none" && it[0].k != "seacells" && it[0].int_cov != "none" }

    sc_binned_groups = COMPUTE_CELL_GROUPS(
        sc_binned_in.map { info, sc_pheno, covs, anno, plink -> tuple(info, covs) }
    )

    sc_binned_full = sc_binned_in
        .combine(sc_binned_groups, by: 0)
        .map { info, sc_pheno, covs, anno, plink, cg ->
            tuple(info, sc_pheno, covs, anno, plink, cg)
        }

    sc_quasar_input_full = sc_no_groups_full
        .mix(sc_seacells_full)
        .mix(sc_binned_full)
        .filter({!(it[0].chr in ["chr1"])})

    sc_quasar = RUN_QUASAR_SC(sc_quasar_input_full)

    sc_quasar_keyed = sc_quasar.map { info, region, variant, time ->
      tuple([cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac], info, region, variant, time)
    }
    sc_quasar = sc_quasar_keyed
        .combine(gene_properties, by: 0)
        .map { base_info, info, region, variant, time, gp_file ->
            tuple(info, region, variant, time, gp_file)
        }

    sc_quasar = COMPUTE_QUASAR_POWER_SC(sc_quasar)
    sc_quasar = COMPUTE_CONVERGENCE_SC(sc_quasar)

    //sc_logcounts = EXTRACT_SC_LOGCOUNTS(
    //    sc_input.filter { it[0].cell_frac == 1.0d && it[0].indiv_frac == 1.0d }
    //)
    //genotype_dosages = EXTRACT_GENOTYPES(all_bed)

    //sc_logcounts_for_figures = sc_logcounts
    //    .map { info, logcounts -> tuple(info.cell_type, logcounts) }
    //sc_int_figures_input = join_sc_pseudo_out
    //    .map { info, covs -> tuple(info.cell_type, info, covs) }
    //    .combine(sc_logcounts_for_figures, by: 0)
    //    .combine(genotype_dosages)
    //    .map { cell_type, info, covs, logcounts, geno ->
    //        tuple(info, covs, logcounts, geno)
    //    }

    pb_quasar_file = Channel
        .of("model\tcell_type\tchr\tcell_frac\tindiv_frac\tint_cov\tpb_type\tregion_file\tvariant_file\ttime_file\tpower_file\tconv_file\tgene_prop_file")
        .concat(pb_quasar
            .map { info, region_file, variant_file, time_file, gene_prop_file, power_file, conv_file  ->
                "${info.model}\t${info.cell_type}\t${info.chr}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${pb_type}\t${region_file}\t${variant_file}\t${time_file}\t${power_file}\t${conv_file}\t${gene_prop_file}"
            }
        )
        .collectFile(name: 'pb_quasar_file', newLine: true, sort: false)
    
    sc_quasar_file = Channel
        .of("cell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tpower_file\tconv_file\tgene_prop_file")
        .concat(sc_quasar
            .map { info, region_file, variant_file, time_file, gene_prop_file, power_file, conv_file ->
                "${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${power_file}\t${conv_file}\t${gene_prop_file}"
            }
        )
        .collectFile(name: 'sc_quasar_file', newLine: true, sort: false)

    seacells_file = Channel
        .of("cell_type\tgroups_file\tinfo_file")
        .concat(seacells_groups
            .combine(seacells_out.sizes, by: 0)
            .map { info, cg, sizes ->
                "${info.cell_type}\t${cg}\t${sizes}"
            }
        )
        .collectFile(
            name: 'seacells_file',
            newLine: true,
            sort: false,
            storeDir: "${projectDir}/data/intermediate"
        )

    // Clumping analysis.
    //clumped_pb_quasar = CLUMP_VARIANTS_PB(pb_quasar.filter{it[0].int_cov == "none"}, all_bed)
    //clumped_pb_quasar_file = Channel
    //    .of("model\tcell_type\tchr\tcell_frac\tindiv_frac\tint_cov\tpb_type\tregion_file\tvariant_file\ttime_file\tpower_file\tconv_file\tgene_prop_file\tclumped_file")
    //    .concat(clumped_pb_quasar
    //        .map { info, region_file, variant_file, time_file, gene_prop_file, power_file, conv_file, clumped_file ->
    //            "${info.model}\t${info.cell_type}\t${info.chr}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${pb_type}\t${region_file}\t${variant_file}\t${time_file}\t${power_file}\t${conv_file}\t${gene_prop_file}\t${clumped_file}"
    //        }
    //    )
    //    .collectFile(
    //        name: 'clumped_pb_quasar_file',
    //        newLine: true,
    //        sort: false,
    //        storeDir: "${projectDir}/data/intermediate"
    //    )

    //clumped_sc_quasar = CLUMP_VARIANTS_SC(sc_quasar.filter{it[0].int_cov == "none"}, all_bed)
    //clumped_sc_quasar_file = Channel
    //    .of("cell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tpower_file\tconv_file\tgene_prop_file\tclumped_file")
    //    .concat(clumped_sc_quasar
    //        .map { info, region_file, variant_file, time_file, gene_prop_file, power_file, conv_file, clumped_file ->
    //            "${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${power_file}\t${conv_file}\t${gene_prop_file}\t${clumped_file}"
    //        }
    //    )
    //    .collectFile(
    //        name: 'clumped_sc_quasar_file',
    //        newLine: true,
    //        sort: false,
    //        storeDir: "${projectDir}/data/intermediate"
    //    )

    // GWAS analysis.

    sc_quasar_gwas_input = sc_quasar_input
        .filter({it[0].cov_spec == "bulk_pca"})
        .filter({it[0].indiv_frac == 1.0d && it[0].cell_frac == 1.0d})
        .filter({it[0].cell_type == "B IN"})
        .combine(Channel.from(1..22))
        .map { info, sc_pheno, covs, anno, bed, pheno_chr ->
          tuple(info + [pheno_chr: pheno_chr as int], sc_pheno, covs, anno, bed)
        }
        .filter({it[0].pheno_chr == 1.0d})

    sc_quasar_gwas = RUN_QUASAR_SC_GWAS(sc_quasar_gwas_input)

    pb_quasar_gwas_input = pb_quasar_input
        .filter({it[0].indiv_frac == 1.0d && it[0].cell_frac == 1.0d})
        .filter({it[0].cell_type == "B IN"})
        .filter({it[0].int_cov == "none"})
        .filter({it[0].model == "nb_glm"})
        .combine(Channel.from(1..22))
        .map { info, pheno, covs, bed, grm, pheno_chr ->
          tuple(info + [pheno_chr: pheno_chr as int], pheno, covs, bed)
        }
        .filter({it[0].pheno_chr == 1.0d})

    pb_quasar_gwas = RUN_QUASAR_PB_GWAS(pb_quasar_gwas_input)

    sc_quasar_gwas_file = Channel
        .of("cell_type\tgeno_chr\tpheno_chr\tsig_variant_file\ttime_file\tn_variants")
        .concat(sc_quasar_gwas
            .map { info, sig_variant_file, time_file, n_variants_file ->
                "${info.cell_type}\t${info.chr}\tchr${info.pheno_chr}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(name: 'sc_quasar_gwas_file', newLine: true, sort: false)

    pb_quasar_gwas_file = Channel
        .of("cell_type\tmodel\tgeno_chr\tpheno_chr\tsig_variant_file\ttime_file\tn_variants")
        .concat(pb_quasar_gwas
            .map { info, sig_variant_file, time_file, n_variants_file ->
                "${info.cell_type}\t${info.model}\t${info.chr}\tchr${info.pheno_chr}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(name: 'pb_quasar_gwas_file', newLine: true, sort: false)

    csaqtl_quasar_file = Channel
        .of("cell_type\tgeno_chr\tsig_variant_file\ttime_file\tn_variants")
        .concat(csaqtl_quasar
            .map { info, sig_variant_file, time_file, n_variants_file ->
                "${info.cell_type}\t${info.chr}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(name: 'csaqtl_quasar_file',newLine: true, sort: false)

    // Genotype permutation analysis.
    rep_bed_files = bed_files
      //.combine(channel.of(1..10))
      .combine(channel.of(1..2))
    permute_bed_files = PERMUTE_BED(rep_bed_files)

    perm_sc_quasar_input = sc_quasar_input_full
      .filter( {it[0].cell_type in ["Plasma", "B IN", "CD4 NC", "B_all", "CD4_T_all", "CD8_T_all"]})
      .filter( {it[0].cov_spec == "bulk_pca" })
      .filter( {it[0].int_cov in ["none", "pseudotime"] + starcat_int_covs})
      .filter( {it[0].k in ["none", 5]})
      .map { info, sc_pheno, covs, anno, bed, cg ->
        tuple(info.chr, info, sc_pheno, covs, anno, bed, cg)
       }
      .combine(permute_bed_files, by: 0)
      .map { chr, info, sc_pheno, covs, anno, bed, cg, perm_bed ->
        tuple(info, sc_pheno, covs, anno, perm_bed, cg)
      }
    
    perm_sc_quasar = RUN_QUASAR_SC_PERM(perm_sc_quasar_input)
    perm_sc_quasar_filt = perm_sc_quasar
      .map { info, region_file, variant_file, time_file ->
        tuple(info.chr, info, region_file, variant_file, time_file)
      }
      .combine(pruned_snps, by: 0)
      .map { chr, info, region_file, variant_file, time_file, prune_in ->
        tuple(info, region_file, variant_file, time_file, prune_in)
      }
      | FILTER_VARIANTS_SC

    perm_sc_quasar_filt = perm_sc_quasar_filt
       .map { info, region_file, variant_file, time_file ->
         tuple([cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac], info, region_file, variant_file, time_file)
       }
       .combine(gene_properties, by: 0)
       .map { base_info, info, region_file, variant_file, time_file, prop_file ->
         tuple(info, region_file, variant_file, time_file, prop_file)
       }

    perm_pb_quasar_input = pb_quasar_input
      .filter( {it[0].cell_type in ["Plasma", "B IN", "CD4 NC", "B_all"]})
      .filter( {it[0].int_cov in ["none", "pseudotime"]})
      .filter( {it[0].model in ["nb_glm", "lmm", "lm"]})
      .map { info, pheno, covs, bed, grm ->
        tuple(info.chr, info, pheno, covs, bed, grm)
       }
      .combine(permute_bed_files, by: 0)
      .map { chr, info, pheno, covs, bed, grm, perm_bed ->
        tuple(info, pheno, covs, perm_bed, grm)
      }

    perm_pb_quasar = RUN_QUASAR_PB_PERM(perm_pb_quasar_input)

    perm_pb_quasar_filt = perm_pb_quasar
      .map { info, region_file, variant_file, time_file ->
        tuple(info.chr, info, region_file, variant_file, time_file)
      }
      .combine(pruned_snps, by: 0)
      .map { chr, info, region_file, variant_file, time_file, prune_in ->
        tuple(info, region_file, variant_file, time_file, prune_in)
      }
      | FILTER_VARIANTS_PB

    perm_pb_quasar_filt = perm_pb_quasar_filt
       .map { info, region_file, variant_file, time_file ->
         tuple([cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac], info, region_file, variant_file, time_file)
       }
       .combine(gene_properties, by: 0)
       .map { base_info, info, region_file, variant_file, time_file, prop_file ->
         tuple(info, region_file, variant_file, time_file, prop_file)
       }

    // GWAS permutations.
    
    gwas_rep_bed_files = bed_files
        .combine(Channel.value(1))
    gwas_permute_bed_files = PERMUTE_BED_GWAS(gwas_rep_bed_files)

    pb_quasar_gwas_perm_input = pb_quasar_input
        .filter({it[0].indiv_frac == 1.0d && it[0].cell_frac == 1.0d})
        .filter({it[0].cell_type == "B IN"})
        .filter({it[0].int_cov == "none"})
        .filter({it[0].model == "nb_glm"})
        .combine(Channel.from(1..22))
        .map { info, pheno, covs, bed, grm, pheno_chr ->
          tuple(info + [pheno_chr: pheno_chr as int], pheno, covs, bed)
        }
        .map { info, pheno, covs, bed ->
          tuple(info.chr, info, pheno, covs, bed)
        }
        .combine(gwas_permute_bed_files, by: 0)
        .map({ it -> [it[1], it[2], it[3], it[5]]})
        .filter({it[0].pheno_chr== 1.0d})

    pb_quasar_gwas_perm = RUN_QUASAR_PB_GWAS_PERM(pb_quasar_gwas_perm_input)

    sc_quasar_gwas_perm_input = sc_quasar_input
        .filter({it[0].cov_spec == "bulk_pca"})
        .filter({it[0].indiv_frac == 1.0d && it[0].cell_frac == 1.0d})
        .filter({it[0].cell_type == "B IN"})
        .combine(Channel.from(1..22))
        .map { info, sc_pheno, covs, anno, bed, pheno_chr ->
         tuple(info + [pheno_chr: pheno_chr as int], sc_pheno, covs, anno, bed)
        }
        .map { info, sc_pheno, covs, anno, bed ->
          tuple(info.chr, info, sc_pheno, covs, anno, bed)
        }
        .combine(gwas_permute_bed_files, by: 0)
        .map({ it -> [it[1], it[2], it[3], it[4], it[6]]})
        .filter({it[0].pheno_chr== 1.0d})

    sc_quasar_gwas_perm = RUN_QUASAR_SC_GWAS_PERM(sc_quasar_gwas_perm_input)

    sc_quasar_gwas_perm_file = Channel
        .of("cell_type\tgeno_chr\tpheno_chr\tsig_variant_file\ttime_file\tn_variants")
        .concat(sc_quasar_gwas_perm
            .map { info, sig_variant_file, time_file, n_variants_file ->
                "${info.cell_type}\t${info.chr}\tchr${info.pheno_chr}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(
            name: 'sc_quasar_gwas_perm_file',
            newLine: true,
            sort: false,
            storeDir: "${projectDir}/data/intermediate"
        )

    pb_quasar_gwas_perm_file = Channel
        .of("cell_type\tmodel\tgeno_chr\tpheno_chr\tsig_variant_file\ttime_file\tn_variants")
        .concat(pb_quasar_gwas_perm
            .map { info, sig_variant_file, time_file, n_variants_file ->
                "${info.cell_type}\t${info.model}\t${info.chr}\tchr${info.pheno_chr}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(
            name: 'pb_quasar_gwas_perm_file',
            newLine: true,
            sort: false,
            storeDir: "${projectDir}/data/intermediate"
        )

    // Permutations under an interaction null.

    rep_pb_cov_files = pb_covs
      .filter({it[0].cell_frac == 1.0d && it[0].indiv_frac == 1.0d})
      .combine(int_cov.filter({it != "none"}))
      .combine(channel.of(1..10))
      .map { info, covs, int_cov, ind ->
        tuple(int_cov, info.cell_type, covs, ind)
      }

    permute_pb_covs = PERMUTE_COV(rep_pb_cov_files)

    perm_int_pb_quasar_input = pb_quasar_input
       .filter( {it[0].cell_frac == 1.0d && it[0].indiv_frac == 1.0d} )
       .filter( {it[0].cell_type in ["Plasma"]} )
       .filter( {it[0].model in ["nb_glm", "lmm", "lm"]} )
       .map { info, pheno, covs, bed, grm ->
        tuple(info.int_cov, info.cell_type, info, pheno, covs, bed, grm)
       }
       .combine(permute_pb_covs, by: [0, 1])
       .map { int_cov, cell_type, info, pheno, covs, bed, grm, perm_cov ->
        tuple(info, pheno, perm_cov, bed, grm)
       }

    perm_int_pb_quasar = RUN_QUASAR_PB_PERM_INT(perm_int_pb_quasar_input)

    perm_int_pb_quasar_filt = perm_int_pb_quasar
      .map { info, region_file, variant_file, time_file ->
        tuple(info.chr, info, region_file, variant_file, time_file)
      }
      .combine(pruned_snps, by: 0)
      .map { chr, info, region_file, variant_file, time_file, prune_in ->
       tuple(info, region_file, variant_file, time_file, prune_in)
      }
      | FILTER_VARIANTS_PB_INT

    perm_int_pb_quasar_filt = perm_int_pb_quasar_filt
       .map { info, region_file, variant_file, time_file ->
         tuple([cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac], info, region_file, variant_file, time_file)
       }
       .combine(gene_properties, by: 0)
       .map { base_info, info, region_file, variant_file, time_file, prop_file ->
         tuple(info, region_file, variant_file, time_file, prop_file)
       }

    perm_int_pb_quasar_file = Channel
        .of("cell_type\tmodel\tchr\tcell_frac\tindiv_frac\tint_cov\tpb_type\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_int_pb_quasar_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.cell_type}\t${info.model}\t${info.chr}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${pb_type}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(
            name: 'pb_quasar_perm_int_file',
            newLine: true,
            sort: false,
            storeDir: "${projectDir}/data/intermediate"
        )

    rep_sc_cov_files = join_sc_pseudo_out
      .combine(channel.of(1..2))
      .map { info, covs, ind ->
        tuple("pseudotime", info.cell_type, covs, ind)
      }

    permute_sc_covs = PERMUTE_SC_COV(rep_sc_cov_files)

    perm_int_sc_quasar_input = sc_quasar_input_full
      .filter { it[0].cell_type == "B_all" && it[0].int_cov == "pseudotime" && it[0].k == "none" }
      .map { info, sc_pheno, covs, anno, bed, cg ->
        tuple(info.int_cov, info.cell_type, info, sc_pheno, covs, anno, bed, cg)
      }
      .combine(permute_sc_covs, by: [0, 1])
      .map { int_cov, cell_type, info, sc_pheno, covs, anno, bed, cg, perm_cov ->
        tuple(info + [interaction_cov: "${info.int_cov}_perm"], sc_pheno, perm_cov, anno, bed, cg)
      }

    perm_int_sc_quasar = RUN_QUASAR_SC_PERM_INT(perm_int_sc_quasar_input)

    perm_int_sc_quasar_filt = perm_int_sc_quasar
      .map { info, region_file, variant_file, time_file ->
        tuple(info.chr, info, region_file, variant_file, time_file)
      }
      .combine(pruned_snps, by: 0)
      .map { chr, info, region_file, variant_file, time_file, prune_in ->
        tuple(info, region_file, variant_file, time_file, prune_in)
      }
      | FILTER_VARIANTS_SC_INT

    perm_int_sc_quasar_filt = perm_int_sc_quasar_filt
       .map { info, region_file, variant_file, time_file ->
         tuple([cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac], info, region_file, variant_file, time_file)
       }
       .combine(gene_properties, by: 0)
       .map { base_info, info, region_file, variant_file, time_file, prop_file ->
         tuple(info, region_file, variant_file, time_file, prop_file)
       }

    sc_quasar_perm_int_file = Channel
        .of("cell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_int_sc_quasar_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(
            name: 'sc_quasar_perm_int_file',
            newLine: true,
            sort: false,
            storeDir: "${projectDir}/data/intermediate"
        )

    rep_sc_cov_files_within = join_sc_pseudo_out
      .combine(channel.of(1..2))
      .map { info, covs, ind ->
        tuple("pseudotime", info.cell_type, covs, ind)
      }

    permute_sc_covs_within = PERMUTE_SC_COV_WITHIN(rep_sc_cov_files_within)

    perm_int_sc_quasar_within_input = sc_quasar_input_full
      .filter { it[0].cell_type == "B_all" && it[0].int_cov == "pseudotime" && it[0].k == "none" }
      .map { info, sc_pheno, covs, anno, bed, cg ->
        tuple(info.int_cov, info.cell_type, info, sc_pheno, covs, anno, bed, cg)
      }
      .combine(permute_sc_covs_within, by: [0, 1])
      .map { int_cov, cell_type, info, sc_pheno, covs, anno, bed, cg, perm_cov ->
        tuple(info + [interaction_cov: "${info.int_cov}_perm"], sc_pheno, perm_cov, anno, bed, cg)
      }

    perm_int_sc_quasar_within = RUN_QUASAR_SC_PERM_WITHIN_INT(perm_int_sc_quasar_within_input)

    perm_int_sc_quasar_within_filt = perm_int_sc_quasar_within
      .map { info, region_file, variant_file, time_file ->
        tuple(info.chr, info, region_file, variant_file, time_file)
      }
      .combine(pruned_snps, by: 0)
      .map { chr, info, region_file, variant_file, time_file, prune_in ->
        tuple(info, region_file, variant_file, time_file, prune_in)
      }
      | FILTER_VARIANTS_SC_INT_WITHIN

    perm_int_sc_quasar_within_filt = perm_int_sc_quasar_within_filt
       .map { info, region_file, variant_file, time_file ->
         tuple([cell_type: info.cell_type, cell_frac: info.cell_frac, indiv_frac: info.indiv_frac], info, region_file, variant_file, time_file)
       }
       .combine(gene_properties, by: 0)
       .map { base_info, info, region_file, variant_file, time_file, prop_file ->
         tuple(info, region_file, variant_file, time_file, prop_file)
       }

    sc_quasar_perm_int_within_file = Channel
        .of("cell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_int_sc_quasar_within_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(
            name: 'sc_quasar_perm_int_within_file',
            newLine: true,
            sort: false,
            storeDir: "${projectDir}/data/intermediate"
        )

    perm_sc_quasar_file = Channel
        .of("cell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_sc_quasar_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(
            name: 'sc_quasar_perm_file',
            newLine: true,
            sort: false,
            storeDir: "${projectDir}/data/intermediate"
        )

    perm_pb_quasar_file = Channel
        .of("cell_type\tmodel\tchr\tcell_frac\tindiv_frac\tint_cov\tpb_type\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_pb_quasar_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.cell_type}\t${info.model}\t${info.chr}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${pb_type}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(
            name: 'pb_quasar_perm_file',
            newLine: true,
            sort: false,
            storeDir: "${projectDir}/data/intermediate"
        )

    // Analysis. 

    //PLOT_CLUMPED(clumped_pb_quasar_file, clumped_sc_quasar_file)
    //PLOT_POWER(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    //PLOT_CONVERGENCE(pb_quasar_file, sc_quasar_file)
    //PLOT_PERM(perm_sc_quasar_file)
    //PLOT_PERM_PB(perm_pb_quasar_file)
    //PLOT_PERM_GLOBAL(perm_pb_quasar_file)
    //PLOT_PERM_GWAS(sc_quasar_gwas_perm_file, pb_quasar_gwas_perm_file)
    //PLOT_PERM_INT(perm_int_pb_quasar_file)
    //PLOT_INT_OUTPUT(pb_quasar_file)
    //PLOT_SC_INT_OUTPUT(sc_quasar_file)
    //PLOT_GROUPED_SC_OUTPUT(sc_quasar_file)
    //PLOT_PERM_SC_GROUPED(perm_sc_quasar_file)
    //PLOT_PERM_SC_INT_GLOBAL(perm_sc_quasar_file)
    //PLOT_PERM_SC_INT(sc_quasar_perm_int_file)
    //PLOT_SC_INT_FIGURES(sc_int_figures_input)
    //PLOT_METACELL_OUTPUT(sc_quasar_file, seacells_file)
    //PLOT_CSAQTL_OUTPUT(csaqtl_quasar_file)
    //PLOT_PERM_SC_INT_WITHIN(sc_quasar_perm_int_within_file)
    //PLOT_TIME(pb_quasar_file, sc_quasar_file, saigeqtl_file, saige_qtl_supp_tables)
    //PLOT_CONCORDANCE(pb_quasar_file, sc_quasar_file, saigeqtl_file)
    //PLOT_GENE_PROPERTIES(gene_properties_file)
    //PLOT_GWAS_OUTPUT(sc_quasar_gwas_file, pb_quasar_gwas_file)
    //PLOT_PVALUE_SCATTER(pb_quasar_file, sc_quasar_file)

    // Create example data for quasar.
    //example_data_input = sc_quasar_input
    //    .filter({it[0].cell_type == "B IN"})
    //    .filter({it[0].indiv_frac == 1.0d && it[0].cell_frac == 1.0d})
    //    .filter({it[0].chr == "chr22"})
    //CREATE_EXAMPLE_DATA(example_data_input)
}
