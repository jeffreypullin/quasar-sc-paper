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
}
