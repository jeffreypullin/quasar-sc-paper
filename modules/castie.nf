
process COLLATE_CASTIE_INPUT {
  label "mega_mem"

  input: tuple val(info), val(sc_counts), val(covs), val(bed), val(anno)
  output: tuple val(info), path("${info.dataset}-${info.cell_type.replaceAll(/\s+/, '-')}-${info.chr}-${info.data_type}-castie-input.tsv"), val(bed), val(anno)

  script:
  def safe = info.cell_type.replaceAll(/\s+/, '-')
  """
  collate-castie-input.R \
    "$info.cell_type" \
    "$info.chr" \
    "$sc_counts" \
    "$covs" \
    "$anno" \
    "$safe" \
    "$info.int_cov"
  mv "${safe}-${info.chr}-castie-input.tsv" "${info.dataset}-${safe}-${info.chr}-${info.data_type}-castie-input.tsv"
  """
}

process RUN_CASTIE {
    label "castie"

    input: tuple val(info), val(input), val(chr_bed), val(anno), val(gene_list), val(subset_bed)
    output: tuple val(info), path("${info.dataset}-rationalised-${info.cell_type}-${info.chr}-${info.data_type}-*-castie-files.tsv")

    script:
    def chr_prefix = "${chr_bed.getParent().toString() + '/' + chr_bed.getSimpleName()}"
    def subset_prefix = "${subset_bed.getParent().toString() + '/' + subset_bed.getSimpleName()}"
    """
    run-castie.py \
        --cell-type "$info.cell_type" \
        --chrom  "$info.chr" \
        --gene-list "$gene_list" \
        --input "$input" \
        --subset-prefix "$subset_prefix" \
        --chr-prefix "$chr_prefix" \
        --anno "$anno" \
        --data-type "$info.data_type" \
        --int-cov "$info.int_cov" \
        --castie-dir "$projectDir/CASTIE" \
        --n-threads ${task.cpus}
    rationalise-castie-files.R *-castie-files.tsv
    for f in rationalised-${info.cell_type}-${info.chr}-${info.data_type}-*-castie-files.tsv; do mv "\$f" "${info.dataset}-\$f"; done
    """
}

process CONCAT_CASTIE_TSVS {

    input: val(tsvs)
    output: path("castie-files.tsv")

    script:
    """
    collate-tsvs.R $tsvs
    mv saigeqtl-files.tsv castie-files.tsv
    """
}
