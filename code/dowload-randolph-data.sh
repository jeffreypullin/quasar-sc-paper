#!/bin/bash

cd /home/jp2045/quasar-sc-paper/data/randolph-data
wget -c -O inputs.tar.gz https://zenodo.org/records/4273999/files/inputs.tar.gz

tar -xzf inputs.tar.gz \
  inputs/1_calculate_pseudobulk/B_cluster_singlets.rds \
  inputs/1_calculate_pseudobulk/CD4_T_cluster_singlets.rds \
  inputs/1_calculate_pseudobulk/CD8_T_cluster_singlets.rds \
  inputs/1_calculate_pseudobulk/monocytes_cluster_singlets.rds \
  inputs/1_calculate_pseudobulk/infected_monocytes_cluster_singlets.rds \
  inputs/1_calculate_pseudobulk/NK_cluster_singlets.rds \
  inputs/1_calculate_pseudobulk/NK_high_response_cluster_singlets.rds \
  inputs/2_calculate_residuals_and_DE_analyses/individual_meta_data_for_GE_with_scaledCovars_with_CTC.txt \
  inputs/3_eQTL_mapping/genotypes.txt \
  inputs/3_eQTL_mapping/SNP_positions.txt \
  inputs/3_eQTL_mapping/GRCh38.92_gene_positions.txt \
  inputs/3_eQTL_mapping/individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt