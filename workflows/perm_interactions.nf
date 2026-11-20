include {
    PERMUTE_COV ;
    PERMUTE_COV as PERMUTE_SC_COV ;
    PERMUTE_SC_COV_WITHIN ;
    PERMUTE_SC_COV_LOGCOUNT_BINS
} from '../modules/permute'
include {
    COMPUTE_CELL_GROUPS as COMPUTE_PERM_CELL_GROUPS_GLOBAL ;
    COMPUTE_CELL_GROUPS as COMPUTE_PERM_CELL_GROUPS_WITHIN
} from '../modules/cellgroups'
include {
    RUN_QUASAR_PB as RUN_QUASAR_PB_PERM_INT ;
    RUN_QUASAR_SC as RUN_QUASAR_SC_PERM_INT ;
    RUN_QUASAR_SC as RUN_QUASAR_SC_PERM_WITHIN_INT ;
    RUN_QUASAR_SC as RUN_QUASAR_SC_PERM_LOGCOUNT_INT
} from '../modules/quasar'
include {
    FILTER_VARIANTS as FILTER_VARIANTS_PB_INT ;
    FILTER_VARIANTS as FILTER_VARIANTS_SC_INT ;
    FILTER_VARIANTS as FILTER_VARIANTS_SC_INT_WITHIN ;
    FILTER_VARIANTS as FILTER_VARIANTS_SC_INT_LOGCOUNT
} from '../modules/postprocess'

// Interaction-null permutations (PB and SC). Extracted from main.nf so the
// entry script stays under Groovy's 65535-character string limit.
workflow PERM_INTERACTIONS {
    take:
    pb_covs
    int_cov
    pb_quasar_input
    pruned_snps
    gene_properties
    sc_covs
    sc_quasar_input_full
    sc_counts

    main:
    def is_starcat_int_cov = { int_cov_name ->
        int_cov_name != null && int_cov_name.toString().startsWith("starcat_")
    }
    def int_cov_for_cov_spec = [
        "bulk_pca+pseudotime": "pseudotime",
        "bulk_pca+starcat_all": "starcat_all",
        "bulk_pca+starcat_CD4_Naive": "starcat_CD4_Naive",
    ]

    // Permutations under an interaction null.
    rep_pb_cov_files = pb_covs
      .filter({it[0].is_full_data})
      .combine(int_cov.filter({it != "none"}))
      .combine(channel.of(1..5))
      // TEMP: skip interaction perms while only running genotype perms.
      .filter { false }
      .map { info, covs, int_cov, ind ->
        tuple(info.dataset, int_cov, info.cell_type, covs, ind)
      }

    permute_pb_covs = PERMUTE_COV(rep_pb_cov_files)

    perm_int_pb_quasar_input = pb_quasar_input
       .filter( {it[0].is_full_data} )
       .filter( {it[0].cell_type in ["Plasma", "B_IN"]} )
       .filter( {it[0].model in ["nb_glm", "nb_glmm", "lmm", "lm", "p_glmm"]} )
       // TEMP: skip interaction perms while only running genotype perms.
       .filter { false }
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

    rep_sc_cov_files = sc_covs
      .filter { info, cov_spec, cov ->
          int_cov_for_cov_spec.containsKey(cov_spec) && !is_starcat_int_cov(int_cov_for_cov_spec[cov_spec])
      }
      .map { info, cov_spec, cov -> tuple(info.dataset, int_cov_for_cov_spec[cov_spec], info.cell_type, cov) }
      .combine(channel.of(1..10))
      // TEMP: GWAS-only run.
      .filter { false }
      .map { dataset, int_cov, cell_type, covs, ind -> tuple(dataset, int_cov, cell_type, covs, ind) }

    permute_sc_covs = PERMUTE_SC_COV(rep_sc_cov_files)

    // Join key for one grouped perm replicate. chr and the perm filename keep
    // replicates and chromosomes from crossing.
    def perm_group_key = { info, perm_cov ->
        [info.dataset, info.chr, info.cell_type, info.model, info.data_type,
         info.cov_spec, info.int_cov, info.k,
         info.cell_frac, info.indiv_frac, info.count_frac, info.min_n_cells,
         info.ablate ?: "none", new File(perm_cov.toString()).name].join("|")
    }

    perm_int_sc_joined = sc_quasar_input_full
      .filter { it[0].data_type != "sct_counts" }
      .filter { it[0].model == "p_glmm_sc" }
      .filter { it[0].k in ["none", 5] }
      .filter { it[0].int_cov != "none" && !is_starcat_int_cov(it[0].int_cov) }
      .map { info, sc_pheno, covs, anno, bed, cg ->
        tuple(info.dataset, info.int_cov, info.cell_type, info, sc_pheno, covs, anno, bed, cg)
      }
      .combine(permute_sc_covs, by: [0, 1, 2])

    // Continuous interaction: swap in the permuted covariate and test G x *_perm.
    // K = 5: QuASAR does not take --interaction, so the null is new groups ranked
    // on the permuted column. Keep the bulk_pca covariates already on the tuple.
    perm_int_sc_quasar_none = perm_int_sc_joined
      .filter { it[3].k == "none" }
      .map { dataset, int_cov, cell_type, info, sc_pheno, covs, anno, bed, cg, perm_cov ->
        tuple(info + [interaction_cov: "${int_cov}_perm"], sc_pheno, perm_cov, anno, bed, cg)
      }

    perm_int_sc_grouped_src = perm_int_sc_joined
      .filter { it[3].k == 5 && it[3].is_full_data }

    perm_int_sc_grouped_groups = COMPUTE_PERM_CELL_GROUPS_GLOBAL(
      perm_int_sc_grouped_src.map { dataset, int_cov, cell_type, info, sc_pheno, covs, anno, bed, cg, perm_cov ->
        tuple(info + [group_cov: "${int_cov}_perm", perm_key: perm_group_key(info, perm_cov)], perm_cov)
      }
    )

    perm_int_sc_quasar_grouped = perm_int_sc_grouped_src
      .map { dataset, int_cov, cell_type, info, sc_pheno, covs, anno, bed, cg, perm_cov ->
        tuple(perm_group_key(info, perm_cov), info, sc_pheno, covs, anno, bed)
      }
      .combine(
        perm_int_sc_grouped_groups.map { ginfo, groups -> tuple(ginfo.perm_key, groups) },
        by: 0
      )
      .map { key, info, sc_pheno, covs, anno, bed, groups ->
        tuple(info, sc_pheno, covs, anno, bed, groups)
      }

    perm_int_sc_quasar_input = perm_int_sc_quasar_none.mix(perm_int_sc_quasar_grouped)

    perm_int_sc_quasar = RUN_QUASAR_SC_PERM_INT(perm_int_sc_quasar_input)
    perm_int_sc_quasar_filt = FILTER_VARIANTS_SC_INT(Utils.combineWithPrunedSnps(perm_int_sc_quasar, pruned_snps))
    perm_int_sc_quasar_filt = Utils.attachGeneProperties(perm_int_sc_quasar_filt, gene_properties)

    rep_sc_cov_files_within = sc_covs
      .filter { info, cov_spec, cov ->
          int_cov_for_cov_spec.containsKey(cov_spec) && !is_starcat_int_cov(int_cov_for_cov_spec[cov_spec])
      }
      .map { info, cov_spec, cov -> tuple(info.dataset, int_cov_for_cov_spec[cov_spec], info.cell_type, cov) }
      .combine(channel.of(1..10))
      // TEMP: GWAS-only run.
      .filter { false }
      .map { dataset, int_cov, cell_type, covs, ind -> tuple(dataset, int_cov, cell_type, covs, ind) }

    permute_sc_covs_within = PERMUTE_SC_COV_WITHIN(rep_sc_cov_files_within)

    perm_int_sc_within_joined = sc_quasar_input_full
      .filter { it[0].data_type != "sct_counts" }
      .filter { it[0].model == "p_glmm_sc" }
      .filter { it[0].k in ["none", 5] }
      .filter { it[0].int_cov != "none" && !is_starcat_int_cov(it[0].int_cov) }
      .map { info, sc_pheno, covs, anno, bed, cg ->
        tuple(info.dataset, info.int_cov, info.cell_type, info, sc_pheno, covs, anno, bed, cg)
      }
      .combine(permute_sc_covs_within, by: [0, 1, 2])

    perm_int_sc_within_none = perm_int_sc_within_joined
      .filter { it[3].k == "none" }
      .map { dataset, int_cov, cell_type, info, sc_pheno, covs, anno, bed, cg, perm_cov ->
        tuple(info + [interaction_cov: "${int_cov}_perm"], sc_pheno, perm_cov, anno, bed, cg)
      }

    perm_int_sc_within_grouped_src = perm_int_sc_within_joined
      .filter { it[3].k == 5 && it[3].is_full_data }

    perm_int_sc_within_grouped_groups = COMPUTE_PERM_CELL_GROUPS_WITHIN(
      perm_int_sc_within_grouped_src.map { dataset, int_cov, cell_type, info, sc_pheno, covs, anno, bed, cg, perm_cov ->
        tuple(info + [group_cov: "${int_cov}_perm", perm_key: perm_group_key(info, perm_cov)], perm_cov)
      }
    )

    perm_int_sc_within_grouped = perm_int_sc_within_grouped_src
      .map { dataset, int_cov, cell_type, info, sc_pheno, covs, anno, bed, cg, perm_cov ->
        tuple(perm_group_key(info, perm_cov), info, sc_pheno, covs, anno, bed)
      }
      .combine(
        perm_int_sc_within_grouped_groups.map { ginfo, groups -> tuple(ginfo.perm_key, groups) },
        by: 0
      )
      .map { key, info, sc_pheno, covs, anno, bed, groups ->
        tuple(info, sc_pheno, covs, anno, bed, groups)
      }

    perm_int_sc_quasar_within_input = perm_int_sc_within_none.mix(perm_int_sc_within_grouped)

    perm_int_sc_quasar_within = RUN_QUASAR_SC_PERM_WITHIN_INT(perm_int_sc_quasar_within_input)
    perm_int_sc_quasar_within_filt = FILTER_VARIANTS_SC_INT_WITHIN(Utils.combineWithPrunedSnps(perm_int_sc_quasar_within, pruned_snps))
    perm_int_sc_quasar_within_filt = Utils.attachGeneProperties(perm_int_sc_quasar_within_filt, gene_properties)

    // Shuffle the interaction covariate within log-library-size bins so that
    // depth-associated pseudotime structure is preserved.
    rep_sc_cov_files_logcount = sc_covs
      .filter { info, cov_spec, cov ->
          int_cov_for_cov_spec.containsKey(cov_spec) && !is_starcat_int_cov(int_cov_for_cov_spec[cov_spec])
      }
      // TEMP: SC interaction perms only on full-data B_all pseudotime.
      .filter { info, cov_spec, cov -> info.is_full_data && info.cell_type == "B_all" }
      // TEMP: GWAS-only run.
      .filter { false }
      .map { info, cov_spec, cov ->
        tuple(
          [info.dataset, info.cell_type, info.cell_frac, info.indiv_frac, info.min_n_cells, info.count_frac],
          info.dataset,
          int_cov_for_cov_spec[cov_spec],
          info.cell_type,
          cov
        )
      }
      .combine(
        sc_counts.map { type, info, counts ->
          tuple(
            [info.dataset, info.cell_type, info.cell_frac, info.indiv_frac, info.min_n_cells, info.count_frac],
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


    emit:
    perm_int_pb_quasar_filt
    perm_int_sc_quasar_filt
    perm_int_sc_quasar_within_filt
    perm_int_sc_quasar_logcount_filt
}
