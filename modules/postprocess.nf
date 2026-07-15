process FILTER_VARIANTS {
    //label "nano"

    input: tuple val(info), val(region_file), val(variant_file), val(time_file), val(prune_in)
    output: tuple val(info), val(region_file), path("${info.dataset}-${info.chr}-${info.cell_type}-quasar-cis-variant-filt.tsv"), val(time_file)

    script:
    """
    filter-variants.R "${info.chr}" "${info.cell_type}" "$variant_file" "$prune_in"
    mv "${info.chr}-${info.cell_type}-quasar-cis-variant-filt.tsv" "${info.dataset}-${info.chr}-${info.cell_type}-quasar-cis-variant-filt.tsv"
    """
}
