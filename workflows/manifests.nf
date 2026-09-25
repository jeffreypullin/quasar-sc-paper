workflow COLLECT_MANIFESTS {
    take:
    gene_properties
    pb_quasar
    sc_quasar
    sc_offset_quasar
    seacells_groups
    seacells_sizes
    sc_quasar_gwas
    pb_quasar_gwas
    csaqtl_quasar
    pc_gwas
    pc_sc_gwas
    sc_quasar_gwas_perm
    pb_quasar_gwas_perm
    perm_int_pb_quasar_filt
    perm_int_sc_quasar_filt
    perm_int_sc_quasar_within_filt
    perm_int_sc_quasar_logcount_filt
    perm_sc_quasar_filt
    perm_pb_quasar_filt

    main:
    gene_properties_file = Channel
        .of("dataset\tcell_type\tproperties_file")
        .concat(
            gene_properties
                .filter { info, properties_file -> info.is_full_data }
                .map { info, properties_file ->
                    "${info.dataset}\t${info.cell_type}\t${properties_file}"
                }
        )
        .collectFile(name: 'gene_properties_file', newLine: true, sort: false)

    pb_quasar_file = Channel
        .of("dataset\tmodel\tcell_type\tchr\tcell_frac\tindiv_frac\tn_cells_target\tcount_frac\tint_cov\tpb_type\tregion_file\tvariant_file\ttime_file\tpower_file\tconv_file\tgene_prop_file")
        .concat(pb_quasar
            .map { info, region_file, variant_file, time_file, gene_prop_file, power_file, conv_file ->
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
            .combine(seacells_sizes, by: 0)
            .map { info, groups, sizes ->
                "${info.dataset}\t${info.cell_type}\t${groups}\t${sizes}"
            }
        )
        .collectFile(name: 'seacells_file', newLine: true, sort: false)

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
        .collectFile(name: 'csaqtl_quasar_file', newLine: true, sort: false)

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
                "${info.dataset}\t${info.cell_type}\t${info.model}\t${info.chr}\tchr${info.pheno_chr}\t${sig_variant_file}\t${time_file}\t${n_variants_file}"
            }
        )
        .collectFile(name: 'pb_quasar_gwas_perm_file', newLine: true, sort: false)

    perm_int_pb_quasar_file = Channel
        .of("dataset\tcell_type\tmodel\tchr\tcell_frac\tindiv_frac\tint_cov\tpb_type\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_int_pb_quasar_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.dataset}\t${info.cell_type}\t${info.model}\t${info.chr}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.pb_type}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(name: 'pb_quasar_perm_int_file', newLine: true, sort: false)

    sc_quasar_perm_int_file = Channel
        .of("dataset\tmodel\tcell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_int_sc_quasar_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.dataset}\t${info.model}\t${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(name: 'sc_quasar_perm_int_file', newLine: true, sort: false)

    sc_quasar_perm_int_within_file = Channel
        .of("dataset\tmodel\tcell_type\tchr\tcov_spec\tcell_frac\tindiv_frac\tint_cov\tk\tregion_file\tvariant_file\ttime_file\tprop_file")
        .concat(perm_int_sc_quasar_within_filt
            .map { info, region_file, variant_file, time_file, prop_file ->
                "${info.dataset}\t${info.model}\t${info.cell_type}\t${info.chr}\t${info.cov_spec}\t${info.cell_frac}\t${info.indiv_frac}\t${info.int_cov}\t${info.k}\t${region_file}\t${variant_file}\t${time_file}\t${prop_file}"
            }
        )
        .collectFile(name: 'sc_quasar_perm_int_within_file', newLine: true, sort: false)

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

    emit:
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
}
