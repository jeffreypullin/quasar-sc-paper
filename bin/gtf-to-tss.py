#!/usr/bin/env python3

import qtl.io

annotation_gtf="Homo_sapiens.GRCh37.82.genes.gtf"
gtf=qtl.io.gtf_to_tss_bed(annotation_gtf, feature='gene')
gtf.to_csv("Homo_sapiens.GRCh37.82.bed.gz",index=False, sep="\t")
