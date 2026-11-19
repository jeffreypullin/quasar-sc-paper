#!/bin/bash
# Download FinnGen R13 summary statistics (+ tabix indexes) for 10
# fairly well-powered autoimmune traits (>= ~3k cases).
set -euo pipefail

OUTDIR=/home/jp2045/quasar-sc-paper/data/finngen-r13-autoimmune
MANIFEST_URL=https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv
PHENOS='M13_RHEUMA|L12_PSORIASIS|K11_IBD_STRICT|T1D_WIDE|K11_UC_STRICT2|K11_COELIAC|E4_GRAVES_STRICT|M13_ANKYLOSPON|M13_SJOGREN|G6_MS'

mkdir -p "$OUTDIR"
cd "$OUTDIR"

wget -c -O finngen_R13_manifest.tsv "$MANIFEST_URL"

awk -F'\t' -v p="$PHENOS" 'NR == 1 || $1 ~ ("^(" p ")$")' finngen_R13_manifest.tsv \
  > selected_manifest.tsv

tail -n +2 selected_manifest.tsv | cut -f7 | while read -r url; do
  wget -c "$url"
  wget -c "${url}.tbi"
done

echo "Downloaded association + .tbi files to $OUTDIR"
gzip -t finngen_R13_*.gz
echo "gzip integrity OK"
