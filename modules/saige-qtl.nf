
process SAIGE_SUBSET_BED {
    conda "$projectDir/envs/cli.yaml"

    input: val plink_bed
    output: path("subset.bed")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """
    plink2 \
      --bfile ${prefix} \
      --mac 20 \
      --write-snplist \
      --out mac20

    shuf mac20.snplist | head -n 1000 > mac20.random1000.snplist

    plink2 \
      --bfile ${prefix} \
      --extract mac20.random1000.snplist \
      --make-bed \
      --out subset
    """
}

process COLLATE_SAIGEQTL_INPUT {
  label "high_mem"

  input: tuple val(info), val(sc_counts), val(expr_pcs), val(geno_pcs), val(bed), val(anno)
  output: tuple val(info), path("${info.cell_type.replaceAll(/\s+/, '-')}-${info.chr}-saigeqtl-input.tsv"), val(bed), val(anno)

  script:
  def safe = info.cell_type.replaceAll(/\s+/, '-')
  """
  collate-saigeqtl-input.R \
    "$info.cell_type" \
    "$info.chr" \
    "$sc_counts" \
    "$expr_pcs" \
    "$geno_pcs" \
    "$anno" \
    "$safe"
  """
}

process EXTRACT_GENES {

    input: 
        tuple val(info), val(input)
        val n
    output: tuple val(info), path("${info.cell_type}-${info.chr}-genes-*.txt")

    script:
    """
    extract-genes.R "$info.cell_type" "$info.chr" "$input" "$n"
    """
}

process RUN_SAIGEQTL {
    label "saigeqtl"

    input: tuple val(info), val(input), val(chr_bed), val(anno), val(gene_list), val(subset_bed)
    output: tuple val(info), path("rationalised-${info.cell_type}-${info.chr}-*-saigeqtl-files.tsv")

    script: 
    def chr_prefix = "${chr_bed.getParent().toString() + '/' + chr_bed.getSimpleName()}"
    def subset_prefix = "${subset_bed.getParent().toString() + '/' + subset_bed.getSimpleName()}"
    """
    run-saigeqtl.py \
        --cell-type "$info.cell_type" \
        --chrom "$info.chr" \
        --gene-list "$gene_list" \
        --input "$input" \
        --subset-prefix "$subset_prefix" \
        --chr-prefix "$chr_prefix" \
        --anno "$anno"
    rationalise-saigeqtl-files.R *-saigeqtl-files.tsv
    """
}

process COMPUTE_SAIGEQTL_POWER {
    label "micro"

    input: tuple val(info), val(files_tsv)
    output: tuple val(info), path("*")

    script:
    """
    compute-saigeqtl-power.R "$files_tsv"
    """
}

process CONCAT_SAIGEQTL_TSVS {
    
    input: val(tsvs)
    output: path("saigeqtl-files.tsv")

    script:
    """
    collate-tsvs.R $tsvs
    """
}
