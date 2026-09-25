process BUILD_TREMOR_SAMPLE_MAP {
    conda "$projectDir/envs/scanpy.yaml"
    label "tiny"

    input: tuple val(dataset), val(vcf), path(curated_map, stageAs: "curated-sample-map.tsv")
    output: tuple val(dataset), path("${dataset}-sample-map.tsv")

    script:
    """
    build-tremor-sample-map.py \\
        "$curated_map" \\
        "$vcf" \\
        "${dataset}-sample-map.tsv"
    """
}

process HARMONISE_TREMOR_H5AD {
    conda "$projectDir/envs/scanpy.yaml"
    label "mega_mem"

    input: tuple val(dataset), val(raw_sc_data), path(sample_map)
    output: tuple val(dataset), path("${dataset}-harmonised.h5ad")

    script:
    """
    harmonise-tremor-h5ad.py $raw_sc_data "${dataset}-harmonised.h5ad" $sample_map
    """
}

process CONVERT_TREMOR_VCF_TO_BED {
    conda "$projectDir/envs/cli.yaml"
    label "tiny"

    input: tuple val(dataset), val(chr), val(vcf), path(sample_map)
    output: tuple val(dataset), val(chr),
        path("${dataset}-${chr}.bed"),
        path("${dataset}-${chr}.bim"),
        path("${dataset}-${chr}.fam")

    shell:
    '''
    gzip -dc !{vcf} | awk '/^#CHROM/ {
      print "#FID\tIID\tSEX"
      for (i = 10; i <= NF; i++) print "0\t" $i "\t2"
      exit
    }' > samples.psam

    # oldFID oldIID newFID newIID
    awk -F'\t' 'NR > 1 { print "0\t" $2 "\t0\t" $1 }' !{sample_map} > update-ids.txt

    plink2 --vcf !{vcf} \
        --psam samples.psam \
        --autosome \
        --maf 0.05 \
        --set-all-var-ids '@:#$r-$a' \
        --make-bed \
        --out tmp \
        --const-fid \
        --output-chr 26

    plink2 --bfile tmp \
        --update-ids update-ids.txt \
        --make-bed \
        --out !{dataset}-!{chr}
    '''
}

process DOWNLOAD_TREMOR_TENSORQTL {
    storeDir "${projectDir}/data/tremor-data/tensorqtl"
    label "nano"

    input: tuple val(dataset), val(cell_type), val(zenodo_id)
    output: tuple val(dataset), val(cell_type), path("${cell_type}.cis_qtl_pairs.zip")

    script:
    """
    curl -fL -o "${cell_type}.cis_qtl_pairs.zip" \\
        "https://zenodo.org/api/records/${zenodo_id}/files/${cell_type}.cis_qtl_pairs.zip/content"
    """
}

process EXTRACT_TREMOR_GENE_MAP {
    conda "$projectDir/envs/tremor-tensorqtl.yaml"
    storeDir "${projectDir}/data/tremor-data"
    label "micro"

    input: tuple val(dataset), path(h5ad)
    output: tuple val(dataset), path("${dataset}-gene-map.tsv")

    script:
    """
    extract-tremor-gene-map.py "$h5ad" "${dataset}-gene-map.tsv"
    """
}

process BUILD_TREMOR_VARIANT_MAP {
    conda "$projectDir/envs/tremor-tensorqtl.yaml"
    storeDir "${projectDir}/data/tremor-data"
    label "micro"

    input: tuple val(dataset), path(vcfs)
    output: tuple val(dataset), path("${dataset}-variant-map.tsv")

    script:
    """
    build-tremor-variant-map.py "${dataset}-variant-map.tsv" $vcfs
    """
}

process BUILD_TREMOR_TQTL_COVS {
    conda "$projectDir/envs/scanpy.yaml"
    label "high_mem"

    input:
        tuple val(info), path(pheno), path(geno_pcs), val(h5ad)
    output:
        tuple val(info), path("${info.dataset}-${info.cell_type}-tqtl-covs.tsv")

    script:
    """
    build-tremor-tqtl-covs.py \\
        "$pheno" \\
        "$geno_pcs" \\
        "$h5ad" \\
        "${info.dataset}-${info.cell_type}-tqtl-covs.tsv"
    """
}

process PREPARE_TREMOR_TENSORQTL_COMPARISON {
    conda "$projectDir/envs/tremor-tensorqtl.yaml"
    label "high_mem"

    input:
        tuple val(dataset),
              val(cell_type),
              path(tensorqtl_zip),
              path(gene_map),
              path(variant_map),
              path(annot_bed),
              path(quasar_variants)
    output:
        tuple val(dataset),
              val(cell_type),
              path("${cell_type}-lm-vs-tensorqtl-sample.tsv"),
              path("${cell_type}-lm-vs-tensorqtl-summary.tsv")

    script:
    """
    prepare-tremor-tensorqtl-comparison.py \\
        --cell-type "${cell_type}" \\
        --tensorqtl-zip "${tensorqtl_zip}" \\
        --gene-map "${gene_map}" \\
        --variant-map "${variant_map}" \\
        --annot-bed "${annot_bed}" \\
        --quasar-variants ${quasar_variants} \\
        --sample-out "${cell_type}-lm-vs-tensorqtl-sample.tsv" \\
        --summary-out "${cell_type}-lm-vs-tensorqtl-summary.tsv"
    """
}
