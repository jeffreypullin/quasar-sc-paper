
process EXTRACT_INDIV_IDS {
    conda "$projectDir/envs/scanpy.yaml"

    input: tuple val(dataset), val(raw_sc_data)
    output: tuple val(dataset), path("${dataset}-indiv-ids.txt")

    script: 
    """
    extract-indiv-ids.py $raw_sc_data
    mv indiv-ids.txt "${dataset}-indiv-ids.txt"
    """
}

process FILTER_VCF {
    conda "$projectDir/envs/cli.yaml"

    input: 
        tuple val(dataset), val(chr), val(vcf), val(sample_file)
    output: tuple val(dataset), val(chr), path("${dataset}-filt-${chr}.vcf.gz"), path("${dataset}-filt-${chr}.vcf.gz.tbi")

    shell:
    """
    bcftools view -S $sample_file $vcf > tmp.vcf
    vcftools --vcf tmp.vcf \
        --maf 0.05 \
        --hwe 1e-6 \
        --recode \
        --recode-INFO-all \
        --stdout \
        --stdout | bgzip -c > "${dataset}-filt-${chr}.vcf.gz"
    tabix -p vcf "${dataset}-filt-${chr}.vcf.gz"
    """
}

process CONVERT_VCF_TO_BED {
    conda "$projectDir/envs/cli.yaml"
    label "tiny"

    input: tuple val(dataset), val(chr), val(vcf), val(vcf_tbi)
    output: tuple val(dataset), val(chr),
        path("${dataset}-${chr}.bed"),
        path("${dataset}-${chr}.bim"),
        path("${dataset}-${chr}.fam")

    shell:
    '''
    plink2 --vcf !{vcf} \
        --set-all-var-ids '@:#$r-$a' \
        --make-bed \
        --out !{dataset}-!{chr} \
        --const-fid
    '''
}

process PRUNE_SNPS {
    conda "$projectDir/envs/cli.yaml"

    input: tuple val(dataset), val(chr), val(plink_bed)
    output: tuple val(dataset), val(chr), path("${dataset}-${chr}.prune.in")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """
    plink2 \
      --bfile "$prefix" \
      --indep-pairwise 250 100 0.3 \
      --rm-dup exclude-mismatch \
      --out "${dataset}-${chr}" \
      --const-fid
    """
}

process CONCAT_BED_FILES {
    conda "$projectDir/envs/cli.yaml"
    label "tiny"

    input: tuple val(dataset), val(bed_files)
    output: tuple val(dataset),
        path("${dataset}-all.bed"),
        path("${dataset}-all.bim"),
        path("${dataset}-all.fam")

    script:
    """
    printf "%s\\n" $bed_files | tr -d '[],' | sort -V -t/ -k9 > all_bed_files.txt
    awk -F. '{print \$1".bed", \$1".bim", \$1".fam"}' all_bed_files.txt > all_files.txt
    plink --keep-allele-order --merge-list all_files.txt --make-bed --out ${dataset}-all
    """
}

process EXTRACT_GENOTYPES {
    conda "$projectDir/envs/extract-genotypes.yaml"
    label "micro"

    input: tuple val(dataset), val(plink_bed), path(variants_tsv)
    output: tuple val(dataset), path("${dataset}-genotype-dosages.tsv")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """
    extract-genotypes.py "$prefix" "$variants_tsv"
    mv genotype-dosages.tsv "${dataset}-genotype-dosages.tsv"
    """
}

process COMPUTE_SC_COUNTS{
    conda "$projectDir/envs/scanpy.yaml"
    label "mega_mem"

    input: tuple val(info), val(raw_sc_data), val(gene_properties)
    output: tuple val("counts"), val(info), path("${info.dataset}-${info.cell_type}-sc-pheno.tsv")

    script:
    """
    compute-sc-counts.py "${info.cell_type}" "${info.cell_frac}" "${info.indiv_frac}" "$raw_sc_data" "$gene_properties" "${info.n_cells_target}" "${info.count_frac}"
    mv "${info.cell_type}-sc-pheno.tsv" "${info.dataset}-${info.cell_type}-sc-pheno.tsv"
    """
}

process COMPUTE_SC_SCT_COUNTS {
    conda "$projectDir/envs/sctransform.yaml"
    label "long_mega_mem"

    input: tuple val(info), path(sc_counts), val(raw_sc_data)
    output: tuple val("sct_counts"), val(info), path("${info.dataset}-${info.cell_type}-sct-sc-pheno.tsv")

    script:
    """
    compute-sc-sct-counts.R "$sc_counts" "$raw_sc_data" "${info.cell_type}"
    mv "${info.cell_type}-sct-sc-pheno.tsv" "${info.dataset}-${info.cell_type}-sct-sc-pheno.tsv"
    """
}

process COMPUTE_SC_OFFSETS {
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"
    publishDir "output", pattern: "*-offset-summary.tsv"

    input: tuple val(info), val(sc_pheno)
    output: tuple val(info),
        path("${info.dataset}-${info.cell_type}-offset-percell.tsv"),
        path("${info.dataset}-${info.cell_type}-offset-donorflat.tsv"),
        path("${info.dataset}-${info.cell_type}-offset-constant.tsv"),
        path("${info.dataset}-${info.cell_type}-offset-summary.tsv")

    script:
    """
    compute-sc-offsets.py "${info.cell_type}" "$sc_pheno"
    for spec in percell donorflat constant summary; do
        mv "${info.cell_type}-offset-\$spec.tsv" "${info.dataset}-${info.cell_type}-offset-\$spec.tsv"
    done
    """
}

process COMPUTE_SC_LOGCOUNTS {
    conda "$projectDir/envs/scanpy.yaml"
    label "mega_mem"

    input: tuple val(info), val(raw_sc_data), val(gene_properties)
    output: tuple val("log_counts"), val(info), path("${info.dataset}-${info.cell_type}-sc-logcounts.tsv")

    script:
    """
    compute-sc-logcounts.py "${info.cell_type}" "${info.cell_frac}" "${info.indiv_frac}" "$raw_sc_data" "$gene_properties" "${info.n_cells_target}" "${info.count_frac}"
    mv "${info.cell_type}-sc-logcounts.tsv" "${info.dataset}-${info.cell_type}-sc-logcounts.tsv"
    """
}

process EXTRACT_SC_LOGCOUNTS {
    conda "$projectDir/envs/scanpy.yaml"
    label "mega_mem"

    input: tuple val(info), val(raw_sc_data), val(gene_properties), path(genes_tsv)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-sc-logcounts.tsv")

    script:
    """
    extract-sc-logcounts.py "${info.cell_type}" "${info.cell_frac}" "${info.indiv_frac}" "$raw_sc_data" "$genes_tsv"
    mv "${info.cell_type}-sc-logcounts.tsv" "${info.dataset}-${info.cell_type}-sc-logcounts.tsv"
    """
}

process COMPUTE_CSAQTL_COUNTS {
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"

    input: tuple val(info), val(groups), val(h5ad)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-csaqtl-pheno.tsv")

    script:
    """
    compute-csaqtl-counts.py "${info.cell_type}" "$groups" "$h5ad"
    mv "${info.cell_type}-csaqtl-pheno.tsv" "${info.dataset}-${info.cell_type}-csaqtl-pheno.tsv"
    """
}

process COMPUTE_SC_PC_PHENO {
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"

    input: tuple val(info), val(raw_sc_data)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-sc-pc-pheno.tsv")

    script:
    """
    compute-sc-pc-pheno.py "${info.cell_type}" "${info.cell_frac}" "${info.indiv_frac}" "$raw_sc_data"
    mv "${info.cell_type}-sc-pc-pheno.tsv" "${info.dataset}-${info.cell_type}-sc-pc-pheno.tsv"
    """
}

process COMPUTE_SC_PC_SC_PHENO {
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"

    input: tuple val(info), val(raw_sc_data)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-sc-pc-sc-pheno.tsv")

    script:
    """
    compute-sc-pc-sc-pheno.py "${info.cell_type}" "${info.cell_frac}" "${info.indiv_frac}" "$raw_sc_data"
    mv "${info.cell_type}-sc-pc-sc-pheno.tsv" "${info.dataset}-${info.cell_type}-sc-pc-sc-pheno.tsv"
    """
}

process COMPUTE_PB_COUNTS{
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"

    input: tuple val(info), val(raw_sc_data), val(gene_properties)
    output: tuple val("counts"), val(info), path("${info.dataset}-${info.cell_type}-pb-pheno.tsv")

    script:
    """
    compute-pb-counts.py "${info.cell_type}" "${info.cell_frac}" "${info.indiv_frac}" "$raw_sc_data" "$gene_properties" "${info.n_cells_target}" "${info.count_frac}"
    mv "${info.cell_type}-pb-pheno.tsv" "${info.dataset}-${info.cell_type}-pb-pheno.tsv"
    """
}

process COMPUTE_PB_LOGCOUNTS{
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"

    input: tuple val(info), val(raw_sc_data), val(gene_properties)
    output: tuple val("logcounts"), val(info), path("${info.dataset}-${info.cell_type}-pb-pheno.tsv")

    script:
    """
    compute-pb-logcounts.py "${info.cell_type}" "${info.cell_frac}" "${info.indiv_frac}" "$raw_sc_data" "$gene_properties" "${info.n_cells_target}" "${info.count_frac}"
    mv "${info.cell_type}-pb-pheno.tsv" "${info.dataset}-${info.cell_type}-pb-pheno.tsv"
    """
}

process COMPUTE_PB_TMM_INT{
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"

    input: tuple val(info), val(raw_sc_data), val(gene_properties)
    output: tuple val("tmm_int"), val(info), path("${info.dataset}-${info.cell_type}-pb-pheno.tsv")

    script:
    """
    compute-pb-tmm-int.py "${info.cell_type}" "${info.cell_frac}" "${info.indiv_frac}" "$raw_sc_data" "$gene_properties" "${info.n_cells_target}" "${info.count_frac}"
    mv "${info.cell_type}-pb-pheno.tsv" "${info.dataset}-${info.cell_type}-pb-pheno.tsv"
    """
}

process COMPUTE_CLUSTER_SIZES{
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"
    publishDir "output"
    
    input: tuple val(dataset), val(raw_sc_data)
    output: tuple val(dataset), path("${dataset}-cluster-sizes.tsv"), path("${dataset}-cells-per-indiv.tsv")

    script:
    """
    compute-cluster-sizes.py "$raw_sc_data"
    mv cluster-sizes.tsv "${dataset}-cluster-sizes.tsv"
    mv cells-per-indiv.tsv "${dataset}-cells-per-indiv.tsv"
    """
}

process COMPUTE_PB_EXPR_COVS{
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"

    input: tuple val(info), val(raw_sc_data), val(gene_properties)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-expr-covs.tsv")

    script:
    """
    compute-pb-expr-covs.py "${info.cell_type}" "$raw_sc_data"
    mv "${info.cell_type}-expr-covs.tsv" "${info.dataset}-${info.cell_type}-expr-covs.tsv"
    """
}

process COMPUTE_SC_EXPR_COVS{
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"

    input: tuple val(info), val(raw_sc_data), val(gene_properties)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-expr-covs.tsv")

    script:
    """
    compute-sc-expr-covs.py "${info.cell_type}" "$raw_sc_data"
    mv "${info.cell_type}-expr-covs.tsv" "${info.dataset}-${info.cell_type}-expr-covs.tsv"
    """
}

process DOWNLOAD_STARCAT_REF {
    conda "$projectDir/envs/starcat.yaml"
    label "long_nano"
    storeDir "$projectDir/data/reference/starcat-cache"

    input: val(url)
    output: path("starcat-cache")

    script:
    """
    mkdir -p starcat-cache
    if [ ! -f starcat-cache/TCAT.V1/TCAT.V1.reference.tsv ]; then
        wget --tries=10 --timeout=120 --waitretry=60 \\
            -O TCAT.V1.tar.gz "$url"
        tar -xzf TCAT.V1.tar.gz -C starcat-cache
        rm TCAT.V1.tar.gz
    fi
    """
}

process COMPUTE_STARCAT_COVS{
    conda "$projectDir/envs/starcat.yaml"
    label "high_mem"
    publishDir "output", pattern: "*.pdf", mode: "copy"

    input: tuple val(info), val(raw_sc_data), val(gene_properties), path(starcat_cache)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-starcat-covs.tsv"), path("${info.dataset}-${info.cell_type}-starcat-score-plots.pdf")

    script:
    """
    compute-starcat-covs.py "${info.cell_type}" "$raw_sc_data" "$starcat_cache"
    mv "${info.cell_type}-starcat-covs.tsv" "${info.dataset}-${info.cell_type}-starcat-covs.tsv"
    mv "${info.cell_type}-starcat-score-plots.pdf" "${info.dataset}-${info.cell_type}-starcat-score-plots.pdf"
    """
}

process COMPUTE_GENOTYPE_PCS {
    conda "$projectDir/envs/cli.yaml"
    label "micro"

    input: tuple val(dataset), val(bed)
    output: tuple val(dataset), path("${dataset}-geno-pcs.txt")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    plink2 --bfile $prefix --indep-pairwise 250 100 0.3 --rm-dup exclude-mismatch --out pruned_variants --threads 2 --const-fid 
    plink2 --bfile $prefix --extract pruned_variants.prune.in --make-bed --out pruned --const-fid 
    plink2 --bfile pruned --pca 10

    cat plink2.eigenvec | \
      sed '1s/IID/sample_id/' | \
      sed '1s/PC/geno_pc/g' | \
      cut -f2- > ${dataset}-geno-pcs.txt
    """
}

process CREATE_ANNOT_BED {
    conda "$projectDir/envs/annotation.yaml"

    input: tuple val(dataset), val(url)
    output: tuple val(dataset), path("${dataset}-annot.bed")

    script:
    """
    wget -O annotation.gtf.gz $url

    collapse-annotation.py \
        annotation.gtf.gz \
        genes.gtf \
        --collapse_only

    gtf-to-tss.py genes.gtf annot.bed.gz

    zcat annot.bed.gz | \
        awk -F'\\t' -v OFS='\\t' '\$1 ~ /(^[1-9]\$)|(^1[0-9]\$)|(^2[012]\$)/ {print \$1,\$2,\$3,\$4}' > \
        ${dataset}-annot.bed
    sed -i "1i #chr\tstart\tend\tphenotype_id" ${dataset}-annot.bed
    """
}

process CREATE_GRM {
    conda "$projectDir/envs/cli.yaml"
    label "micro"
    
    input: tuple val(dataset), val(bed)
    output: tuple val(dataset), path("${dataset}-grm.tsv")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    plink2 --bfile $prefix --indep-pairwise 250 50 0.2 --out onek1k_pruning_info --threads 2
    plink2 --bfile $prefix --extract onek1k_pruning_info.prune.in --make-king square --out king_ibd_out --threads 2
    process-grm.R onek1k
    mv grm.tsv "${dataset}-grm.tsv"
    """
}

process PREPARE_SLINGSHOT_ADATA {
    conda "$projectDir/envs/prepare-slingshot.yaml"
    label "long_mega_mem"

    input: tuple val(info), val(raw_sc_data)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-slingshot-input.tsv.gz")

    script:
    """
    prepare-slingshot-adata.R "${info.cell_type}" "$raw_sc_data"
    mv "${info.cell_type}-slingshot-input.tsv.gz" "${info.dataset}-${info.cell_type}-slingshot-input.tsv.gz"
    """
}

process PREPARE_SEACELLS_ADATA {
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"

    input: tuple val(info), val(raw_sc_data)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-seacells-input.h5ad")

    script:
    """
    prepare-seacells-adata.py "${info.cell_type}" "$raw_sc_data"
    mv "${info.cell_type}-seacells-input.h5ad" "${info.dataset}-${info.cell_type}-seacells-input.h5ad"
    """
}

process RUN_SLINGSHOT {
    conda "$projectDir/envs/slingshot.yaml"
    label "high_mem"
    publishDir path: "output", pattern: "*.pdf", mode: "copy"

    input: tuple val(info), path(slingshot_input)
    output:
        tuple val(info), path("${info.dataset}-${info.cell_type}-pseudotime.tsv"),
        path("${info.dataset}-${info.cell_type}-slingshot.pdf"),
        path("${info.dataset}-${info.cell_type}-slingshot-celltype.pdf")

    script:
    """
    run-slingshot.R "$slingshot_input" "${info.cell_type}"
    mv "${info.cell_type}-pseudotime.tsv" "${info.dataset}-${info.cell_type}-pseudotime.tsv"
    mv "${info.cell_type}-slingshot.pdf" "${info.dataset}-${info.cell_type}-slingshot.pdf"
    mv "${info.cell_type}-slingshot-celltype.pdf" "${info.dataset}-${info.cell_type}-slingshot-celltype.pdf"
    """
}

process JOIN_SC_INT_COVS {
    label "micro"

    input: tuple val(info), path(expr_covs), path(int_cov)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-cov-data-with-int-cov.tsv")

    script:
    """
    join-sc-int-covs.R "${info.cell_type}" "$expr_covs" "$int_cov" 
    mv "${info.cell_type}-cov-data-with-int-cov.tsv" "${info.dataset}-${info.cell_type}-cov-data-with-int-cov.tsv"
    """
}

process JOIN_PB_INT_COVS {
    label "micro"

    input: tuple val(info), path(pb_covs), path(sc_expr_covs), path(int_cov)
    output: tuple val(info), path("${info.dataset}-${info.cell_type}-pb-cov-data-with-int-cov.tsv")

    script:
    """
    join-pb-int-covs.R "${info.cell_type}" "$pb_covs" "$sc_expr_covs" "$int_cov"
    mv "${info.cell_type}-pb-cov-data-with-int-cov.tsv" "${info.dataset}-${info.cell_type}-pb-cov-data-with-int-cov.tsv"
    """
}
