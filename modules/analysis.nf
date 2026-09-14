
process PLOT_POWER {
    publishDir "output"

    input: 
        val pb_quasar_file
        val sc_quasar_file
        val saigeqtl_file
    output: tuple path("power-method-comparison-plot.pdf"), 
        path("power-frac-plot.pdf"),
        path("power-sc-non-zero-frac-plot.pdf")

    script:
    """
    plot-power.R $pb_quasar_file $sc_quasar_file $saigeqtl_file
    """
}

process PLOT_POWER_N_CELLS_FILTER {
    publishDir "output"

    input:
        val pb_quasar_file
        val sc_quasar_file
    output:
        path "power-n-cells-filter-plot.pdf"

    script:
    """
    plot-power-n-cells-filter.R $pb_quasar_file $sc_quasar_file
    """
}

process PLOT_POWER_OFFSET {
    publishDir "output"
    label "high_mem"

    input:
        val sc_offset_quasar_file
    output:
        tuple path("power-offset-plot.pdf"),
              path("power-offset-scatter-plot.pdf"),
              path("power-offset-summary.tsv")

    script:
    """
    plot-power-offset.R $sc_offset_quasar_file
    """
}

process PLOT_CELLS_PER_INDIV {
    publishDir "output"

    input:
        path cells_per_indiv_file
    output:
        path "cells-per-indiv-plot.pdf"

    script:
    """
    plot-cells-per-indiv.R $cells_per_indiv_file
    """
}

process PLOT_POWER_COVS {
    publishDir "output"

    input:
        val sc_quasar_file
    output:
        tuple path("power-cov-spec-plot.pdf"),
              path("power-cov-spec-gene-plot.pdf"),
              path("power-int-cov-plot.pdf")

    script:
    """
    plot-power-covs.R $sc_quasar_file
    """
}

process COMPUTE_CONVERGENCE {
    label "nano"
    
    input: tuple val(info), val(region_file), val(variant_file), val(time_file), val(gene_properties_files), val(power_file)
    output: tuple val(info), val(region_file), val(variant_file), val(time_file), val(gene_properties_files), val(power_file),
        path("${info.dataset}-${info.chr}-${info.cell_type}-problem-variants.tsv")

    script:
    """
    compute-convergence.R "${info.dataset}" "${info.chr}" "${info.cell_type}" "${info.int_cov}" "$variant_file" "$gene_properties_files"
    """    
}

process CLUMP_VARIANTS {
    label "long_nano"

    input:
        tuple val(info), val(region_file), val(variant_file), val(time_file), val(gene_prop_file), val(power_file), val(conv_file), val(all_bed)
    output: tuple val(info), val(region_file), val(variant_file), val(time_file), val(gene_prop_file), val(power_file), val(conv_file), 
        path("${info.dataset}-${info.chr}-${info.cell_type}-${info.int_cov}-harmonised-clumps.tsv")

    script:
    def prefix = all_bed.getParent().toString() + '/' + all_bed.getSimpleName()
    """
    prepare-clump-inputs.R "$variant_file"

    for assoc in assoc-*.tsv; do
      fid="\${assoc#assoc-}"
      fid="\${fid%.tsv}"
      plink2 \
        --bfile "${prefix}" \
        --clump "\$assoc" \
        --clump-p1 5e-6 \
        --clump-p2 5e-2 \
        --clump-r2 0.5 \
        --out "${info.chr}-${info.cell_type}-${info.int_cov}-clump-\${fid}"
    done

    harmonise-clump-files.R "${info.chr}" "${info.cell_type}" "${info.int_cov}"
    mv "${info.chr}-${info.cell_type}-${info.int_cov}-harmonised-clumps.tsv" "${info.dataset}-${info.chr}-${info.cell_type}-${info.int_cov}-harmonised-clumps.tsv"
    """
}

process PLOT_CONVERGENCE {
    publishDir "output"
    label "high_mem"

    input: 
        val pb_quasar_file
        val sc_quasar_file
        val saigeqtl_file
    output: path("convergence-plot.pdf")

    script:
    """
    plot-convergence.R $pb_quasar_file $sc_quasar_file $saigeqtl_file
    """
}

process PLOT_PERM {
    label "high_mem"
    publishDir "output"

    input: val perm_sc_quasar_file
    output: path("perm-overall-plot.pdf")

    script:
    """
    plot-perm.R $perm_sc_quasar_file
    """
}

process PLOT_PERM_CELL_FRAC {
    label "high_mem"
    publishDir "output"

    input:
        val perm_sc_quasar_file
        val perm_pb_quasar_file
    output:
        tuple path("perm-cell-frac-lm-plot.pdf"),
              path("perm-cell-frac-nb_glm-plot.pdf"),
              path("perm-cell-frac-p_glmm_sc-plot.pdf"),
              path("perm-cell-frac-lmm_sc-plot.pdf")

    script:
    """
    plot-perm-cell-frac.R $perm_sc_quasar_file $perm_pb_quasar_file
    """
}

process PLOT_PERM_COUNT_FRAC {
    label "high_mem"
    publishDir "output"

    input:
        val perm_sc_quasar_file
        val perm_pb_quasar_file
    output:
        tuple path("perm-count-frac-lm-plot.pdf"),
              path("perm-count-frac-nb_glm-plot.pdf"),
              path("perm-count-frac-p_glmm_sc-plot.pdf"),
              path("perm-count-frac-lmm_sc-plot.pdf")

    script:
    """
    plot-perm-count-frac.R $perm_sc_quasar_file $perm_pb_quasar_file
    """
}

process PLOT_PERM_PB {
    publishDir "output"
    label "high_mem"

    input: val perm_pb_quasar_file
    output: tuple path("perm-pb-overall-plot.pdf"),
        path("perm-pb-cell-frac-plot.pdf")

     script:
    """
    plot-perm-pb.R $perm_pb_quasar_file
    """
}

process PLOT_PERM_GLOBAL {
    publishDir "output"
    label "high_mem"

    input: val perm_pb_quasar_file
    output: tuple path("perm-global-pb-overall-int-plot.pdf"), path("perm-global-pb-overall-int-main-plot.pdf")

    script:
    """
    plot-perm-global.R $perm_pb_quasar_file
    """
}

process PLOT_PERM_INT {
    publishDir "output"
    label "high_mem"

    input: val perm_pb_quasar_file
    output: tuple path("perm-int-pb-overall-int-plot.pdf"), path("perm-int-pb-overall-int-main-plot.pdf")

    script:
    """
    plot-perm-int.R $perm_pb_quasar_file
    """
}

process PLOT_PERM_GWAS {
    publishDir "output"
    label "high_mem"

    input: 
      val perm_gwas_sc_quasar_file
      val perm_gwas_pb_quasar_file
    output: path("perm-gwas-overall-plot.pdf")
    
    script:
    """
    plot-perm-gwas.R $perm_gwas_sc_quasar_file $perm_gwas_pb_quasar_file
    """

}

process PLOT_INT_OUTPUT {
    publishDir "output"
    label "high_mem"

    input: val perm_pb_quasar_file
    output: path("int-res-egenes.tsv")

    script:
    """
    plot-int-output.R $perm_pb_quasar_file
    """
}

process PLOT_PERM_SC_INT_GLOBAL {
    publishDir "output"
    label "high_mem"

    input: val perm_sc_quasar_file
    output: tuple path("perm-sc-int-global-plot.pdf"), path("perm-sc-int-global-main-plot.pdf")

    script:
    """
    plot-perm-sc-int-global.R $perm_sc_quasar_file
    """
}

process PLOT_PERM_SC_INT {
    publishDir "output"
    label "high_mem"

    input: val perm_sc_quasar_file
    output: tuple path("perm-sc-int-plot.pdf"), path("perm-sc-int-main-plot.pdf")

    script:
    """
    plot-perm-sc-int.R $perm_sc_quasar_file
    """
}

process PLOT_PERM_SC_INT_WITHIN {
    publishDir "output"
    label "high_mem"

    input: val perm_sc_quasar_file
    output: tuple path("perm-sc-int-within-plot.pdf"), path("perm-sc-int-within-main-plot.pdf")

    script:
    """
    plot-perm-sc-int-within.R $perm_sc_quasar_file
    """
}

process PLOT_PERM_SC_INT_LOGCOUNT_BINS {
    publishDir "output"
    label "high_mem"

    input: val perm_sc_quasar_file
    output: tuple path("perm-sc-int-logcount-bins-plot.pdf"),
                 path("perm-sc-int-logcount-bins-main-plot.pdf")

    script:
    """
    plot-perm-sc-int-logcount-bins.R $perm_sc_quasar_file
    """
}

process PLOT_PERM_SC_GROUPED {
    publishDir "output"
    label "high_mem"

    input: val perm_sc_quasar_file
    output: path("perm-sc-grouped-*-plot.pdf")

    script:
    """
    plot-perm-sc-grouped.R $perm_sc_quasar_file
    """
}

process PLOT_GROUPED_SC_OUTPUT {
    label "high_mem"
    publishDir "output"

    input: val sc_quasar_file
    output: tuple path("grouped-sc-output-*-linear-egenes-plot.pdf"),
                 path("grouped-sc-output-*-acat-egenes-plot.pdf"),
                 path("grouped-sc-output-het-vs-linear-plot.pdf"),
                 path("grouped-sc-output-gws-leads.tsv")

    script:
    """
    plot-grouped-sc-output.R $sc_quasar_file
    """
}

process PLOT_GROUPED_VS_INT {
    label "high_mem"
    publishDir "output"

    input:
        val sc_quasar_file
        val castie_file
    output: tuple path("grouped-vs-int-scatter-plot.pdf"),
                 path("grouped-vs-int-upset-plot.pdf"),
                 path("grouped-vs-int-upset-genes.tsv")

    script:
    """
    plot-grouped-vs-int.R $sc_quasar_file $castie_file
    """
}

process PLOT_GROUPED_INT_TIME {
    publishDir "output"

    input: val sc_quasar_file
    output: path("grouped-int-time-plot.pdf")

    script:
    """
    plot-grouped-int-time.R $sc_quasar_file
    """
}

process PLOT_SC_PGLMM_VS_LMM {
    label "high_mem"
    publishDir "output"

    input: val sc_quasar_file
    output: path("sc-pglmm-vs-lmm-plot.pdf")

    script:
    """
    plot-sc-pglmm-vs-lmm.R $sc_quasar_file
    """
}

process PLOT_SC_INT_FIGURES {
    publishDir "output"

    input: tuple val(info), path(sc_covs), path(sc_logcounts), path(genotype_dosages)
    output: path("sc-int-figures.pdf")

    script:
    """
    plot-sc-int-figures.R "$sc_covs" "$sc_logcounts" "$genotype_dosages"
    """
}

process PLOT_METACELL_OUTPUT {
    label "high_mem"
    publishDir "output"

    input:
        val sc_quasar_file
        val seacells_file
    output: path("metacell-output-sizes-plot.pdf")

    script:
    """
    plot-metacell-output.R $sc_quasar_file $seacells_file
    """
}

process PLOT_SC_INT_OUTPUT {
    label "high_mem"
    publishDir "output"

    input:
        val sc_quasar_file
        val castie_file
    output: path("sc-int-eqtl-hits.tsv")

    script:
    """
    plot-sc-int-output.R $sc_quasar_file $castie_file
    """
}

process PLOT_QUASAR_SC_INT_OUTPUT {
    label "high_mem"
    publishDir "output"

    input: val sc_quasar_file
    output: path("quasar-sc-int-eqtl-hits.tsv")

    script:
    """
    plot-sc-int-output.R $sc_quasar_file
    mv sc-int-eqtl-hits.tsv quasar-sc-int-eqtl-hits.tsv
    """
}

process PLOT_GWAS_OUTPUT {
    publishDir "output"

    input: 
      val sc_quasar_gwas_file
      val pb_quasar_gwas_file
    output: path("gwas-sig-counts-plot.pdf")

    script:
    """
    plot-gwas-output.R $sc_quasar_gwas_file $pb_quasar_gwas_file
    """
}

process PLOT_CSAQTL_OUTPUT {
    publishDir "output"

    input: val csaqtl_quasar_file
    output: path("csaqtl-output-plot.pdf")

    script:
    """
    plot-csaqtl-output.R $csaqtl_quasar_file
    """
}

process PLOT_PC_GWAS_OUTPUT {
    label "high_mem"
    publishDir "output"

    input:
        val pc_gwas_file
        val pc_sc_gwas_file
    output:
        tuple path("pc-gwas-B_IN-plot.pdf"),
        path("pc-gwas-CD4_NC-plot.pdf"),
        path("pc-gwas-lmm_sc-B_IN-plot.pdf"),
        path("pc-gwas-lmm_sc-CD4_NC-plot.pdf"),
        path("pc-gwas-gws-leads.tsv")

    script:
    """
    plot-pc-gwas-output.R $pc_gwas_file $pc_sc_gwas_file
    """
}

process PLOT_PC_GWAS_COMPARISON {
    publishDir "output"

    input:
        val pc_gwas_file
        val pc_sc_gwas_file
    output:
        tuple path("time-pc-gwas-plot.pdf"),
              path("power-pc-gwas-plot.pdf")

    script:
    """
    plot-pc-gwas-comparison.R $pc_gwas_file $pc_sc_gwas_file
    """
}

process PLOT_PVALUE_SCATTER {
    publishDir "output"
    label "high_mem"

    input:
      val pb_quasar_file
      val sc_quasar_file

    output: path("pvalue-scatter.pdf")

    script:
    """
    plot-pvalue-scatter.R $pb_quasar_file $sc_quasar_file
    """
}

process PLOT_ZSCORE_SCATTER {
    publishDir "output"
    label "high_mem"

    input:
      val pb_quasar_file
      val sc_quasar_file
      val saigeqtl_file

    output:
      path("zscore-scatter-B_IN.pdf")
      path("zscore-scatter-saigeqtl-B_IN.pdf")

    script:
    """
    plot-zscore-scatter.R $pb_quasar_file $sc_quasar_file $saigeqtl_file
    """
}

process PLOT_METHOD_SCATTER {
    publishDir "output"
    label "high_mem"

    input:
        val pb_quasar_file
        val sc_quasar_file
    output:
        path("method-scatter-plot.pdf")

    script:
    """
    plot-method-scatter.R $pb_quasar_file $sc_quasar_file
    """
}

process PLOT_MAIN_VS_INT {
    publishDir "output"
    label "high_mem"

    input: val sc_quasar_file

    output: path("main-vs-int-*-plot.pdf")

    script:
    """
    plot-main-vs-int.R $sc_quasar_file
    """
}

process PLOT_INT_TIME {
    publishDir "output"
    label "high_mem"

    input:
        val sc_quasar_file
        val castie_file
    output: path("time-int-plot.pdf")

    script:
    """
    plot-int-time.R $sc_quasar_file $castie_file
    """
}

process PLOT_QUASAR_INT_TIME {
    publishDir "output"
    label "high_mem"

    input: val sc_quasar_file
    output: path("time-int-plot.pdf")

    script:
    """
    plot-int-time.R $sc_quasar_file
    """
}

process PLOT_TIME {
    publishDir "output"

    input: 
        val pb_quasar_file
        val sc_quasar_file
        val saigeqtl_file
    output: tuple path("time-quasar-plot.pdf"), 
        path("time-method-comparison-plot.pdf"),
        path("time-frac-plot.pdf")

    script:
    """
    plot-time.R $pb_quasar_file $sc_quasar_file $saigeqtl_file
    """ 
}

process PLOT_CONCORDANCE {
    publishDir "output"

    input: 
        val pb_quasar_file
        val sc_quasar_file
        val saigeqtl_file
    output: tuple path("concordance-plot.pdf"),
        path("concordance-lm-p_glmm_sc-plot.pdf")

    script:
    """
    plot-concordance.R $pb_quasar_file $sc_quasar_file $saigeqtl_file
    """
}

process PLOT_UNIQUE_SC_EGENES {
    publishDir "output"
    label "high_mem"

    input:
        val pb_quasar_file
        val sc_quasar_file
    output:
        path("unique-p-glmm-egenes.tsv")

    script:
    """
    plot-unique-sc-egenes.R $pb_quasar_file $sc_quasar_file
    """
}

process PLOT_UNIQUE_SC_EGENE_FIGURES {
    publishDir "output"
    label "high_mem"

    input:
        tuple val(info), path(leads_tsv), path(sc_logcounts), path(genotype_dosages)
    output:
        path("unique-sc-egene-figures-${info.cell_type}.pdf")

    script:
    """
    plot-unique-sc-egene-figures.R "$leads_tsv" "$sc_logcounts" "$genotype_dosages" "${info.cell_type}"
    """
}

process PLOT_UNIQUE_PB_EGENES {
    publishDir "output"
    label "high_mem"

    input:
        val pb_quasar_file
        val sc_quasar_file
    output:
        path("unique-lm-egenes.tsv")

    script:
    """
    plot-unique-pb-egenes.R $pb_quasar_file $sc_quasar_file
    """
}

process PLOT_UNIQUE_PB_EGENE_FIGURES {
    publishDir "output"
    label "high_mem"

    input:
        tuple val(info), path(leads_tsv), path(sc_logcounts), path(genotype_dosages)
    output:
        path("unique-pb-egene-figures-${info.cell_type}.pdf")

    script:
    """
    plot-unique-pb-egene-figures.R "$leads_tsv" "$sc_logcounts" "$genotype_dosages" "${info.cell_type}"
    """
}

process PLOT_EGENE_SIG_MODEL {
    publishDir "output"
    label "high_mem"

    input:
        val pb_quasar_file
        val sc_quasar_file
    output:
        tuple path("egene-sig-model-forest-plot.pdf"),
              path("egene-sig-model-logit.tsv")

    script:
    """
    plot-egene-sig-model.R $pb_quasar_file $sc_quasar_file
    """
}

process PLOT_POWER_CELL_FRAC {
    publishDir "output"
    label "high_mem"

    input:
        val pb_quasar_file
        val sc_quasar_file
    output:
        tuple path("power-cell-frac-counts-plot.pdf"),
              path("power-cell-frac-recall-plot.pdf")

    script:
    """
    plot-power-cell-frac.R $pb_quasar_file $sc_quasar_file
    """
}

process PLOT_POWER_COUNT_FRAC {
    publishDir "output"
    label "high_mem"

    input:
        val pb_quasar_file
        val sc_quasar_file
    output:
        tuple path("power-count-frac-counts-plot.pdf"),
              path("power-count-frac-recall-plot.pdf")

    script:
    """
    plot-power-count-frac.R $pb_quasar_file $sc_quasar_file
    """
}

process PLOT_POWER_INDIV_FRAC {
    publishDir "output"
    label "high_mem"

    input:
        val pb_quasar_file
        val sc_quasar_file
    output:
        tuple path("power-indiv-frac-counts-plot.pdf"),
              path("power-indiv-frac-recall-plot.pdf")

    script:
    """
    plot-power-indiv-frac.R $pb_quasar_file $sc_quasar_file
    """
}

process PLOT_POWER_BOTH_FRAC {
    publishDir "output"
    label "high_mem"

    input:
        val pb_quasar_file
        val sc_quasar_file
    output:
        path "power-both-frac-heatmap.pdf"

    script:
    """
    plot-power-both-frac.R $pb_quasar_file $sc_quasar_file
    """
}

process PLOT_POWER_GENE_PROP {
    publishDir "output"
    label "high_mem"

    input:
        val pb_quasar_file
        val sc_quasar_file
    output:
        tuple path("power-gene-prop-sc-mean-plot.pdf"),
              path("power-gene-prop-sc-var-plot.pdf"),
              path("power-gene-prop-sc-cv-plot.pdf"),
              path("power-gene-prop-sc-non-zero-frac-plot.pdf"),
              path("power-gene-prop-pb-mean-plot.pdf"),
              path("power-gene-prop-pb-var-plot.pdf"),
              path("power-gene-prop-pb-cv-plot.pdf"),
              path("power-gene-prop-pb-non-zero-frac-plot.pdf")

    script:
    """
    plot-power-gene-prop.R $pb_quasar_file $sc_quasar_file
    """
}

process PLOT_CLUMPED {
    publishDir "output"

    input:
        val pb_clumped_file
        val sc_clumped_file
    output: tuple path("clumped-n-signals-plot.pdf"),
        path("clumped-multi-signal-plot.pdf"),
        path("clumped-n-independent-eqtls-plot.pdf")

    script:
    """
    plot-clumped.R $pb_clumped_file $sc_clumped_file
    """
}

process COMPUTE_GENE_PROPERTIES {
    label "high_mem"

    input:
        tuple val(info), val(adata)
        val anno_file
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-gene-properties.tsv")

    script:
    """
    compute-gene-properties.py "${info.cell_type}" "${info.cell_frac}" "${info.indiv_frac}" "$adata" "$anno_file" "${info.n_cells_target}" "${info.count_frac}"
    mv "${info.cell_type}-gene-properties.tsv" "${info.dataset}-${info.cell_type}-gene-properties.tsv"
    """
}

process PLOT_GENE_PROPERTIES {
    publishDir "output"

    input: val gene_properties_file
    output: tuple path("sc-non-zero-frac-vs-pb-non-zero-frac-plot.pdf"), 
      path("sc-non-zero-frac-bin-gene-counts-plot.pdf"), 
      path("sc-non-zero-frac-vs-sc-mean-plot.pdf")

    script: 
    """
    plot-gene-properties.R "$gene_properties_file"
    """
}

process CREATE_EXAMPLE_DATA {
    conda "$projectDir/envs/cli.yaml"
    publishDir "output"

    input: tuple val(info), val(pheno_bed), val(covs), val(anno), val(plink_bed)
    output: tuple path("sc-pheno-n100.tsv"), path("example-anno.tsv")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """ 
    awk '{print \$2}' ${prefix}.fam | head -n 100 > first_100_ids.txt
    plink2 --bfile $prefix --keep first_100_ids.txt --make-bed --out chr22-n100
    make-example-data.R first_100_ids.txt "$pheno_bed" "$anno"
    """
}

process RUN_COLOC {
    publishDir "output"

    input: val sc_quasar_file
    output: path("coloc-results.tsv")

    script:
    """
    run-coloc.R $sc_quasar_file
    """
}
