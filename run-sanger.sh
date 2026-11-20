#!/bin/bash
# Launch the pipeline on the Sanger LSF farm (Cytoimmgen by default).

set -euo pipefail

mkdir -p /lustre/scratch124/humgen/projects_v2/healthy_imm_expr/cytoimmgen_quasar/tmp
export NXF_TEMP="/lustre/scratch124/humgen/projects_v2/healthy_imm_expr/cytoimmgen_quasar/tmp"

args=("$@")

# Default to cytoimmgen when no --datasets flag is passed.
has_datasets=0
for a in "${args[@]+"${args[@]}"}"; do
  case "$a" in
    --datasets|--datasets=*) has_datasets=1 ;;
  esac
done
if [[ "$has_datasets" -eq 0 ]]; then
  args+=(--datasets cytoimmgen)
fi

./nextflow -log /dev/null run main.nf -profile sanger -resume "${args[@]}"
