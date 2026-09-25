process CONVERT_RANDOLPH_SEURAT {
    conda "$projectDir/envs/randolph-seurat.yaml"
    label "high_mem"

    input:
        tuple val(dataset),
              val(cluster),
              path(seurat_rds),
              path(metadata)
    output:
        tuple val(dataset),
              val(cluster),
              path("${cluster}-mtx")

    script:
    """
    convert-randolph-seurat.R \\
        "$seurat_rds" \\
        "$metadata" \\
        "$cluster" \\
        "${cluster}-mtx"
    """
}

process ASSEMBLE_RANDOLPH_H5AD {
    conda "$projectDir/envs/scanpy.yaml"
    label "mega_mem"

    input:
        tuple val(dataset),
              val(condition),
              path(mtx_dirs)
    output:
        tuple val("${dataset}_${condition}"),
              path("${dataset}_${condition}.h5ad")

    script:
    def mtx_args = mtx_dirs.collect { "\"${it}\"" }.join(" \\\n        ")
    """
    assemble-randolph-h5ad.py \\
        "$condition" \\
        "${dataset}_${condition}.h5ad" \\
        $mtx_args
    """
}

process CONVERT_RANDOLPH_GENOTYPES {
    conda "$projectDir/envs/cli.yaml"
    label "randolph_genotypes"
    storeDir "${projectDir}/data/randolph-data/plink"

    input:
        tuple val(dataset),
              path(genotypes),
              path(snp_positions)
    output:
        tuple val(dataset),
              path("${dataset}-chr*.bed"),
              path("${dataset}-chr*.bim"),
              path("${dataset}-chr*.fam")

    script:
    """
    convert-randolph-genotypes.py \\
        "$genotypes" \\
        "$snp_positions" \\
        "${dataset}"
    """
}

process CONVERT_RANDOLPH_ANNOT {
    conda "$projectDir/envs/scanpy.yaml"
    label "micro"
    storeDir "${projectDir}/data/randolph-data"

    input:
        tuple val(dataset),
              path(gene_positions)
    output:
        tuple val(dataset),
              path("${dataset}-annot.bed")

    script:
    """
    convert-randolph-annot.py \\
        "$gene_positions" \\
        "${dataset}-annot.bed"
    """
}

process BUILD_RANDOLPH_COVS {
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"

    input:
        tuple val(info),
              path(pheno),
              path(geno_pcs),
              path(metadata),
              path(ctc_metadata)
    output:
        tuple val(info),
              path("${info.dataset}-${info.cell_type}-covs.tsv")

    script:
    """
    build-randolph-covs.py \\
        "$pheno" \\
        "$geno_pcs" \\
        "$metadata" \\
        "$ctc_metadata" \\
        "${info.cell_type}" \\
        "${info.condition}" \\
        "${info.dataset}-${info.cell_type}-covs.tsv"
    """
}

process EXPAND_RANDOLPH_SC_COVS {
    conda "$projectDir/envs/scanpy.yaml"
    label "micro"

    input:
        tuple val(info),
              path(pb_covs),
              path(sc_pheno)
    output:
        tuple val(info),
              path("${info.dataset}-${info.cell_type}-sc-covs.tsv")

    script:
    """
    expand-randolph-sc-covs.py \\
        "$pb_covs" \\
        "$sc_pheno" \\
        "${info.dataset}-${info.cell_type}-sc-covs.tsv"
    """
}
