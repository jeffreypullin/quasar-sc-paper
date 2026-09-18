
process COMBINE_PB_COVS {

    input: 
        tuple val(info), val(expr_covs), val(geno_pcs)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-covs.tsv")

    script:
    """
    combine-pb-covs.R "${info.cell_type}" "$expr_covs" "$geno_pcs"
    mv "${info.cell_type}-covs.tsv" "${info.dataset}-${info.cell_type}-covs.tsv"
    """
}

process COMBINE_SC_COVS {
    label "micro"

    input: 
        tuple val(info), val(pb_expr_covs), val(sc_expr_covs), val(geno_pcs)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-covs.tsv")

    script:
    """
    combine-sc-covs.R "${info.cell_type}" "$pb_expr_covs" "$sc_expr_covs" "$geno_pcs"
    mv "${info.cell_type}-covs.tsv" "${info.dataset}-${info.cell_type}-covs.tsv"
    """
}

process FILTER_COVS {
    input: tuple val(info), val(cov), val(cov_spec)
    output: tuple val(info), val(cov_spec), path("${info.dataset}-${info.cell_type}-${cov_spec}-covs.tsv")
    """
    filter-covs.R "${info.cell_type}" "$cov_spec" "$cov"
    mv "${info.cell_type}-${cov_spec}-covs.tsv" "${info.dataset}-${info.cell_type}-${cov_spec}-covs.tsv"
    """
}

process ANNOTATE_PHENO {

    input: tuple val(pb_type), val(info), val(expr_covs), path(annot_bed)
    output: tuple val(info), val(pb_type), path("${info.dataset}-${info.cell_type}-annot-pheno.tsv")

    script:
    """
    annotate-pheno.R "${info.cell_type}" "$expr_covs" "$annot_bed"
    mv "${info.cell_type}-annot-pheno.tsv" "${info.dataset}-${info.cell_type}-annot-pheno.tsv"
    """
}

process RUN_QUASAR_PB {
    label "pb_quasar"

    input: tuple val(info), val(pheno_bed), val(covs), val(plink_bed), val(grm)
    output: tuple val(info),
        path("${info.dataset}-${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}-quasar-cis-region.txt.gz"), 
        path("${info.dataset}-${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}-quasar-cis-variant.txt.gz"),
        path("${info.dataset}-${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}-time.txt")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    def apl_flag = (info.model in ["nb_glm", "nb_glmm"]) ? "--use-apl" : ""
    def grm_flag = (info.model in ["lmm", "p_glmm", "nb_glmm"]) ? "-g ${grm}" : ""
    def quant_res_flag = (info.model == "nb_glm") ? "--use-quant-res" : ""
    def interaction_cov = info.interaction_cov ?: info.int_cov
    def int_flag = (info.int_cov != "none") ? "-i ${interaction_cov}" : ""
    """
    /usr/bin/time -p -o "${info.dataset}-${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      -p "$prefix" \
      -b "$pheno_bed" \
      -c "$covs" \
      -o "${info.dataset}-${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}" \
      --model "${info.model}" \
      --mode                                                                             cis \
      ${quant_res_flag} \
      ${apl_flag} \
      ${int_flag} \
      ${grm_flag} \
      --verbose
    gzip "${info.dataset}-${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}-quasar-cis-variant.txt"
    gzip "${info.dataset}-${info.model}-${info.chr}-${info.cell_type}-${info.int_cov}-quasar-cis-region.txt"
    """
}


process RUN_QUASAR_SC {
    label "sc_quasar"

    input: tuple val(info), val(pheno_bed), val(covs), val(anno), val(plink_bed), val(cell_groups)
    output: tuple val(info),
        path("${info.dataset}-${info.chr}-${info.cell_type}-${info.model}-${info.data_type}-${info.int_cov}-K${info.k}-quasar-cis-region.txt.gz"), 
        path("${info.dataset}-${info.chr}-${info.cell_type}-${info.model}-${info.data_type}-${info.int_cov}-K${info.k}-quasar-cis-variant.txt.gz"),
        path("${info.dataset}-${info.chr}-${info.cell_type}-${info.model}-${info.data_type}-${info.int_cov}-K${info.k}-time.txt")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    def interaction_cov = info.interaction_cov ?: info.int_cov
    def int_flag = (info.int_cov != "none" && info.k == "none") ? "-i '${interaction_cov}'" : ""
    def cg_flag = (info.k != "none") ? "--cell-groups ${cell_groups}" : ""
    def base = "${info.dataset}-${info.chr}-${info.cell_type}-${info.model}-${info.data_type}-${info.int_cov}-K${info.k}"
    """
    /usr/bin/time -p -o "${base}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      --plink "$prefix" \
      --sc-pheno "$pheno_bed" \
      --anno "$anno" \
      --cov "$covs" \
      --out "${base}" \
      --model                                                                                                                                                                                                                                                                 ${info.model} \
      --mode cis \
      ${int_flag} \
      ${cg_flag} \
      --verbose
    gzip "${base}-quasar-cis-variant.txt"
    gzip "${base}-quasar-cis-region.txt"
    """
}

process RUN_QUASAR_SC_OFFSET {
    label "sc_quasar"

    input: tuple val(info), val(pheno_bed), val(covs), val(offset_file), val(anno), val(plink_bed)
    output: tuple val(info),
        path("${info.dataset}-${info.chr}-${info.cell_type}-${info.model}-${info.offset_spec}-quasar-cis-region.txt.gz"),
        path("${info.dataset}-${info.chr}-${info.cell_type}-${info.model}-${info.offset_spec}-quasar-cis-variant.txt.gz"),
        path("${info.dataset}-${info.chr}-${info.cell_type}-${info.model}-${info.offset_spec}-time.txt")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    def base = "${info.dataset}-${info.chr}-${info.cell_type}-${info.model}-${info.offset_spec}"
    """
    /usr/bin/time -p -o "${base}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      --plink "$prefix" \
      --sc-pheno "$pheno_bed" \
      --anno "$anno" \
      --cov "$covs" \
      --offset-file "$offset_file" \
      --out "${base}" \
      --model ${info.model} \
      --mode   cis \
      --verbose
    gzip "${base}-quasar-cis-variant.txt"
    gzip "${base}-quasar-cis-region.txt"
    """
}

process RUN_QUASAR_SC_GWAS {
    label "sc_quasar_gwas"

    input: tuple val(info), val(pheno_bed), val(covs), val(anno), val(plink_bed)
    output: tuple val(info),
        path("sig-${info.dataset}-pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt"),
        path("${info.dataset}-pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-time.txt"),
        path("${info.dataset}-pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-n-variants.txt")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """
    /usr/bin/time -p -o "${info.dataset}-pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      --plink "$prefix" \
      --sc-pheno "$pheno_bed" \
      --anno "$anno" \
      --cov "$covs" \
      --out "${info.dataset}-pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}" \
      --model     ${info.model} \
      --pheno-chr ${info.pheno_chr} \
      --mode gwas \
      --verbose
    
    awk 'BEGIN {FS=OFS="\t"} NR==1 || (\$10+0) < 5e-6' \
        "${info.dataset}-pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt" >\
        "sig-${info.dataset}-pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt"
    awk 'END { print NR-1 }' \
        "${info.dataset}-pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt" >\
        "${info.dataset}-pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-n-variants.txt"
    rm "${info.dataset}-pheno-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt"
    """
}

process RUN_QUASAR_PB_GWAS {
    label "sc_quasar_gwas"

    input: tuple val(info), val(pheno_bed), val(covs), val(plink_bed)
    output: tuple val(info),
        path("sig-${info.dataset}-pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt"),
        path("${info.dataset}-pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-time.txt"),
        path("${info.dataset}-pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-n-variants.txt")

    script:
    def apl_flag = (info.model == "nb_glm") ? "--use-apl" : ""
    def quant_res_flag = (info.model == "nb_glm") ? "--use-quant-res" : ""
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """
    /usr/bin/time -p -o "${info.dataset}-pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      -p "$prefix" \
      -b "$pheno_bed" \
      -c "$covs" \
      --out "${info.dataset}-pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}" \
      --model "${info.model}" \
      --pheno-chr ${info.pheno_chr} \
      --mode   gwas \
      ${apl_flag} \
      ${quant_res_flag} \
      --verbose
    
    awk 'BEGIN {FS=OFS="\t"} NR==1 || (toupper(\$10)!="NAN" && (\$10+0) < 5e-6)' \
        "${info.dataset}-pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt" >\
        "sig-${info.dataset}-pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt"
    awk 'END { print NR-1 }' \
        "${info.dataset}-pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt" >\
        "${info.dataset}-pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-n-variants.txt"
    rm "${info.dataset}-pheno-${info.model}-chr${info.pheno_chr}-geno-${info.chr}-${info.cell_type}-quasar-gwas-variant.txt"
    """
}

process RUN_QUASAR_GWAS {
    label "sc_quasar_gwas"

    input: tuple val(info), val(pheno_bed), val(covs), val(plink_bed)
    output: tuple val(info),
        path("sig-${info.dataset}-${info.gwas_label}-${info.model}-${info.cell_type}-quasar-gwas-variant.txt"),
        path("${info.dataset}-${info.gwas_label}-${info.model}-${info.cell_type}-time.txt"),
        path("${info.dataset}-${info.gwas_label}-${info.model}-${info.cell_type}-n-variants.txt")

    script:
    def out_prefix = "${info.dataset}-${info.gwas_label}-${info.model}-${info.cell_type}"
    def apl_flags = (info.model == 'nb_glm') ? '--use-apl --use-quant-res' : ''
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """
    /usr/bin/time -p -o "${out_prefix}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      -p "$prefix" \
      -b "$pheno_bed" \
      -c "$covs" \
      --out "${out_prefix}" \
      --model ${info.model} \
      --mode   gwas \
      ${apl_flags} \
      --verbose

    awk 'BEGIN {FS=OFS="\t"} NR==1 || (toupper(\$10)!="NAN" && (\$10+0) < 5e-6)' \
        "${out_prefix}-quasar-gwas-variant.txt" >\
        "sig-${out_prefix}-quasar-gwas-variant.txt"
    awk 'END { print NR-1 }' \
        "${out_prefix}-quasar-gwas-variant.txt" >\
        "${out_prefix}-n-variants.txt"
    rm "${out_prefix}-quasar-gwas-variant.txt"
    """
}

process RUN_QUASAR_SC_PC_GWAS {
    label "sc_quasar_gwas"

    input: tuple val(info), val(pheno), val(covs), val(plink_bed)
    output: tuple val(info),
        path("sig-${info.dataset}-${info.gwas_label}-${info.model}-${info.cell_type}-quasar-gwas-variant.txt"),
        path("${info.dataset}-${info.gwas_label}-${info.model}-${info.cell_type}-time.txt"),
        path("${info.dataset}-${info.gwas_label}-${info.model}-${info.cell_type}-n-variants.txt")

    script:
    def out_prefix = "${info.dataset}-${info.gwas_label}-${info.model}-${info.cell_type}"
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """
    /usr/bin/time -p -o "${out_prefix}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      --plink "$prefix" \
      --sc-pheno "$pheno" \
      --cov "$covs" \
      --out "${out_prefix}" \
      --model ${info.model} \
      --mode   gwas \
      --verbose

    awk 'BEGIN {FS=OFS="\t"} NR==1 || (toupper(\$10)!="NAN" && (\$10+0) < 5e-6)' \
        "${out_prefix}-quasar-gwas-variant.txt" >\
        "sig-${out_prefix}-quasar-gwas-variant.txt"
    awk 'END { print NR-1 }' \
        "${out_prefix}-quasar-gwas-variant.txt" >\
        "${out_prefix}-n-variants.txt"
    rm "${out_prefix}-quasar-gwas-variant.txt"
    """
}

process COMPUTE_QUASAR_POWER {
    label "nano"
    
    input: tuple val(info), val(region_file), val(variant_file), val(time_file), val(gene_properties_files)
    output: tuple val(info), val(region_file), val(variant_file), val(time_file), val(gene_properties_files),
        path("${info.dataset}-${info.chr}-${info.cell_type}-n-sig-variants.tsv")

    script:
    """
    compute-quasar-power.R "${info.chr}" "${info.cell_type}" "${info.int_cov}" "$variant_file" "$gene_properties_files"
    mv "${info.chr}-${info.cell_type}-n-sig-variants.tsv" "${info.dataset}-${info.chr}-${info.cell_type}-n-sig-variants.tsv"
    """
}
