process HARMONISE_CYTOIMMGEN_H5AD {
    conda "$projectDir/envs/scanpy.yaml"
    label "mega_mem"

    input:
        // val keeps large lustre inputs in place (same pattern as Tremor).
        tuple val(dataset),
              val(raw_sc_data),
              val(predictions),
              val(obs_cols)
    output:
        tuple val(dataset), path("${dataset}-harmonised.h5ad")

    script:
    def overrides = obs_cols.collect { k, v -> "${k}=${v}" }.join(" ")
    """
    harmonise-cytoimmgen-h5ad.py \\
        "$raw_sc_data" \\
        "$predictions" \\
        "${dataset}-harmonised.h5ad" \\
        $overrides
    """
}

process SPLIT_CYTOIMMGEN_PLINK {
    conda "$projectDir/envs/cli.yaml"
    label "tiny"

    input:
        // val (not path) so plink reads the lustre prefix in place; staging
        // the .bed alone would leave .bim/.fam behind.
        tuple val(dataset),
              val(chr),
              val(plink_bed),
              path(indiv_ids)
    output:
        tuple val(dataset),
              val(chr),
              path("${dataset}-${chr}.bed"),
              path("${dataset}-${chr}.bim"),
              path("${dataset}-${chr}.fam")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    def chr_num = chr.toString().replaceAll(/^chr/, "")
    """
    awk '{ print "0\\t" \$1 }' "$indiv_ids" > keep.txt

    plink2 --bfile "$prefix" \\
        --chr ${chr_num} \\
        --keep keep.txt \\
        --maf 0.05 \\
        --output-chr 26 \\
        --set-all-var-ids '@:#\$r-\$a' \\
        --make-bed \\
        --out "${dataset}-${chr}" \\
        --const-fid
    """
}

process BUILD_CYTOIMMGEN_SC_COVS {
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"

    input:
        tuple val(info),
              path(sc_pheno),
              path(pb_expr_covs),
              path(geno_pcs),
              path(harmonised_h5ad),
              val(cov_spec),
              val(timepoints)
    output:
        tuple val(info),
              val(cov_spec),
              path("${info.dataset}-${info.cell_type}-${cov_spec}-covs.tsv")

    script:
    def tp = timepoints instanceof List ? timepoints.join(",") : timepoints.toString()
    def out = "${info.dataset}-${info.cell_type}-${cov_spec}-covs.tsv"
    """
    build-cytoimmgen-sc-covs.py \\
        "$sc_pheno" \\
        "$pb_expr_covs" \\
        "$geno_pcs" \\
        "$harmonised_h5ad" \\
        "${info.cell_type}" \\
        "$cov_spec" \\
        "$out" \\
        "$tp"
    """
}
