process COMPUTE_CELL_GROUPS {
    //label "nano"

    input: tuple val(info), val(cov_file)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-${info.int_cov}-K${info.k}-cell-groups.tsv")

    script:
    """
    compute-cell-groups.R "$cov_file" "${info.int_cov}" "${info.k}" "${info.cell_type}"
    mv "${info.cell_type}-${info.int_cov}-K${info.k}-cell-groups.tsv" "${info.dataset}-${info.cell_type}-${info.int_cov}-K${info.k}-cell-groups.tsv"
    """
}

process COMPUTE_SEACELLS {
    conda "$projectDir/envs/seacells.yaml"
    label "mega_mem"

    input: tuple val(info), path(h5ad)
    output:
        tuple val(info), path("${info.dataset}-${info.cell_type}-seacells-cell-groups.tsv"), emit: groups
        tuple val(info), path("${info.dataset}-${info.cell_type}-seacells-info.tsv"),        emit: sizes

    script:
    """
    compute-seacells.py "$h5ad" "${info.cell_type}"
    mv "${info.cell_type}-seacells-cell-groups.tsv" "${info.dataset}-${info.cell_type}-seacells-cell-groups.tsv"
    mv "${info.cell_type}-seacells-info.tsv" "${info.dataset}-${info.cell_type}-seacells-info.tsv"
    """
}
