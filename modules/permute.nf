process PERMUTE_BED {

    input: tuple val(dataset), val(chr), val(plink_bed), val(ind)
    output: tuple val(dataset), val(chr), path("permute-${chr}.bed")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """
    permute-fam.R ${prefix}.fam ${chr} ${ind}
    cp ${prefix}.bim ./"permute-${chr}.bim"
    cp ${prefix}.bed ./"permute-${chr}.bed"
    """
}

process PERMUTE_COV {

    input: tuple val(dataset), val(int_cov), val(cell_type), val(cov), val(ind)
    output: tuple val(dataset), val(int_cov), val(cell_type), path("permute-${cell_type}-${int_cov}-${ind}.tsv")

    script:
    """
    permute-cov.R "$ind" "$cell_type" "$int_cov" "$cov"
    """
}

process PERMUTE_SC_COV_WITHIN {

    input: tuple val(dataset), val(int_cov), val(cell_type), val(cov), val(ind)
    output: tuple val(dataset), val(int_cov), val(cell_type), path("permute-within-${cell_type}-${int_cov}-${ind}.tsv")

    script:
    """
    permute-sc-cov-within.R "$ind" "$cell_type" "$int_cov" "$cov"
    """
}
