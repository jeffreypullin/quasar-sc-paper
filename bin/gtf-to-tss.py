#!/usr/bin/env python3

import sys

import qtl.io

annotation_gtf, out_path = sys.argv[1:3]
gtf = qtl.io.gtf_to_tss_bed(annotation_gtf, feature="gene")
gtf.to_csv(out_path, index=False, sep="\t")
