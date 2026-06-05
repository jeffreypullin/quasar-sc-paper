
process COMBINE_PB_COVS {

    input: 
        tuple val(info), val(expr_covs)
        val geno_pcs
    output: tuple val(info), path("${info.cell_type}-covs.tsv")

    script:
    """
    combine-pb-covs.R "${info.cell_type}" "$expr_covs" "$geno_pcs"
    """
}

process COMBINE_SC_COVS {
    label "micro"

    input: 
        tuple val(info), val(pb_expr_covs), val(sc_expr_covs)
        val geno_pcs
    output: tuple val(info), path("${info.cell_type}-covs.tsv")

    script:
    """
    combine-sc-covs.R "${info.cell_type}" "$pb_expr_covs" "$sc_expr_covs" "$geno_pcs"
    """
}

process FILTER_COVS {
    input: tuple val(info), val(cov), val(cov_spec)
    output: tuple val(info), val(cov_spec), path("${info.cell_type}-${cov_spec}-covs.tsv")
    """
    filter-covs.R "${info.cell_type}" "$cov_spec" "$cov"
    """
}

process ANNOTATE_PHENO {

    input: 
        tuple val(pb_type), val(info), val(expr_covs)
        val annot_bed
    output: tuple val(info), val(pb_type), path("${info.cell_type}-annot-pheno.tsv")

    script:
    """
    annotate-pheno.R "${info.cell_type}" "$expr_covs" "$annot_bed"
    """
}

process RUN_QUASAR_PB {
    label "pb_quasar"

    input: tuple val(info), val(pheno_bed), val(covs), val(plink_bed), val(grm)
    output: tuple val(info),
        path("${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}-quasar-cis-region.txt.gz"), 
        path("${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}-quasar-cis-variant.txt.gz"),
        path("${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}-time.txt")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    def apl_flag = (info.model == "nb_glm") ? "--use-apl" : ""
    def grm_flag = (info.model.contains("lmm")) ? "-g ${grm}" : ""
    def quant_res_flag = (info.model == "nb_glm") ? "--use-quant-res" : ""
    def int_flag = (info.int_cov != "none") ? "-i ${info.int_cov}" : ""
    """
    /usr/bin/time -p -o "${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      -p "$prefix" \
      -b "$pheno_bed" \
      -c "$covs" \
      -o "${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}" \
      --model "${info.model}" \
      --mode                        cis \
      ${quant_res_flag} \
      ${apl_flag} \
      ${int_flag} \
      ${grm_flag} \
      --verbose
    gzip "${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}-quasar-cis-variant.txt"
    gzip "${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}-quasar-cis-region.txt"
    """
}

process COMPUTE_CELL_GROUPS {
    //label "nano"

    input: tuple val(info), val(cov_file)
    output: tuple val(info), path("${info.cell_type}-${info.int_cov}-K${info.k}-cell-groups.tsv")

    script:
    """
    compute-cell-groups.R "$cov_file" "${info.int_cov}" "${info.k}" "${info.cell_type}"
    """
}

process COMPUTE_SEACELLS {
    conda "$projectDir/envs/seacells.yaml"
    label "mega_mem"
    publishDir path: "output", pattern: "*.tsv", mode: "copy"

    input: tuple val(info), path(h5ad)
    output:
        tuple val(info), path("${info.cell_type.replaceAll(' ', '_')}-seacells-cell-groups.tsv"), emit: groups
        tuple val(info), path("${info.cell_type.replaceAll(' ', '_')}-seacells-info.tsv"),        emit: sizes

    script:
    """
    compute-seacells.py "$h5ad" "${info.cell_type}"
    """
}

process RUN_QUASAR_SC {
    label "sc_quasar"

    input: tuple val(info), val(pheno_bed), val(covs), val(anno), val(plink_bed), val(cell_groups)
    output: tuple val(info),
        path("${info.chr}-${info.cell_type}-${info.int_cov}-K${info.k}-quasar-cis-region.txt.gz"), 
        path("${info.chr}-${info.cell_type}-${info.int_cov}-K${info.k}-quasar-cis-variant.txt.gz"),
        path("${info.chr}-${info.cell_type}-${info.int_cov}-K${info.k}-time.txt")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    def interaction_cov = info.interaction_cov ?: info.int_cov
    def int_flag = (info.int_cov != "none" && info.k == "none") ? "-i ${interaction_cov}" : ""
    def cg_flag = (info.k != "none") ? "--cell-groups ${cell_groups}" : ""
    def base = "${info.chr}-${info.cell_type}-${info.int_cov}-K${info.k}"
    """
    /usr/bin/time -p -o "${base}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      --plink "$prefix" \
      --sc-pheno "$pheno_bed" \
      --anno "$anno" \
      --cov "$covs" \
      --out "${base}" \
      --model                                                                                    p_glmm_sc \
      --mode                                        cis \
      ${int_flag} \
      ${cg_flag} \
      --verbose
    gzip "${base}-quasar-cis-variant.txt"
    gzip "${base}-quasar-cis-region.txt"
    """
}

process RUN_QUASAR_SC_GWAS {
    label "sc_quasar_gwas"

    input: tuple val(info), val(pheno_bed), val(covs), val(anno), val(plink_bed)
    output: tuple val(info),
        path("sig-pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt"),
        path("pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-time.txt"),
        path("pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-n-variants.txt")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """
    /usr/bin/time -p -o "pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      --plink "$prefix" \
      --sc-pheno "$pheno_bed" \
      --anno "$anno" \
      --cov "$covs" \
      --out "pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}" \
      --model   p_glmm_sc \
      --pheno-chr ${info.pheno_chr} \
      --mode gwas \
      --verbose
    
    awk 'BEGIN {FS=OFS="\t"} NR==1 || (\$10+0) < 5e-6' \
        "pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt" >\
        "sig-pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt"
    awk 'END { print NR-1 }' \
        "pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt" >\
        "pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-n-variants.txt"
    rm "pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt"
    """
}

process RUN_QUASAR_PB_GWAS {
    label "sc_quasar_gwas"

    input: tuple val(info), val(pheno_bed), val(covs), val(plink_bed)
    output: tuple val(info),
        path("sig-pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt"),
        path("pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-time.txt"),
        path("pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-n-variants.txt")

    script:
    def apl_flag = (info.model == "nb_glm") ? "--use-apl" : ""
    def quant_res_flag = (info.model == "nb_glm") ? "--use-quant-res" : ""
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """
    /usr/bin/time -p -o "pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      -p "$prefix" \
      -b "$pheno_bed" \
      -c "$covs" \
      --out "pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}" \
      --model "${info.model}" \
      --pheno-chr ${info.pheno_chr} \
      --mode gwas \
      ${apl_flag} \
      ${quant_res_flag} \
      --verbose
    
    awk 'BEGIN {FS=OFS="\t"} NR==1 || (toupper(\$10)!="NAN" && (\$10+0) < 5e-6)' \
        "pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt" >\
        "sig-pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt"
    awk 'END { print NR-1 }' \
        "pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt" >\
        "pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-n-variants.txt"
    rm "pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt"
    """
}

process RUN_CSAQTL_QUASAR {
    label "sc_quasar_gwas"

    input: tuple val(info), val(pheno_bed), val(covs), val(plink_bed)
    output: tuple val(info),
        path("sig-csaqtl-nb_glm-${info.cell_type}-quasar-gwas-variant.txt"),
        path("csaqtl-nb_glm-${info.cell_type}-time.txt"),
        path("csaqtl-nb_glm-${info.cell_type}-n-variants.txt")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """
    /usr/bin/time -p -o "csaqtl-nb_glm-${info.cell_type}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      -p "$prefix" \
      -b "$pheno_bed" \
      -c "$covs" \
      --out "csaqtl-nb_glm-${info.cell_type}" \
      --model nb_glm \
      --mode gwas \
      --use-apl \
      --use-quant-res \
      --verbose

    awk 'BEGIN {FS=OFS="\t"} NR==1 || (toupper(\$10)!="NAN" && (\$10+0) < 5e-6)' \
        "csaqtl-nb_glm-${info.cell_type}-quasar-gwas-variant.txt" >\
        "sig-csaqtl-nb_glm-${info.cell_type}-quasar-gwas-variant.txt"
    awk 'END { print NR-1 }' \
        "csaqtl-nb_glm-${info.cell_type}-quasar-gwas-variant.txt" >\
        "csaqtl-nb_glm-${info.cell_type}-n-variants.txt"
    rm "csaqtl-nb_glm-${info.cell_type}-quasar-gwas-variant.txt"
    """
}

process FILTER_VARIANTS {
    //label "nano"

    input: tuple val(info), val(region_file), val(variant_file), val(time_file), val(prune_in)
    output: tuple val(info), val(region_file), path("${info.chr}-${info.cell_type}-quasar-cis-variant-filt.tsv"), val(time_file)

    script:
    """
    filter-variants.R "${info.chr}" "${info.cell_type}" "$variant_file" "$prune_in"
    """
}

process COMPUTE_QUASAR_POWER {
    //label "nano"
    
    input: tuple val(info), val(region_file), val(variant_file), val(time_file), val(gene_properties_files)
    output: tuple val(info), val(region_file), val(variant_file), val(time_file), val(gene_properties_files),
        path("${info.chr}-${info.cell_type}-n-sig-variants.tsv")

    script:
    """
    compute-quasar-power.R "${info.chr}" "${info.cell_type}" "${info.int_cov}" "$variant_file" "$gene_properties_files"
    """
}
