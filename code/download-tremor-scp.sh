#!/bin/bash

cd /home/jp2045/quasar-sc-paper/data/tremor-data
# NB: This download url is time limited.
curl "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP3177&auth_code=lczqJA6T&directory=all&context=study"  -o cfg.txt; curl -K cfg.txt && rm cfg.txt
cp SCP3177/documentation/ET.cerebellum_eQTL_fordownload.h5ad ET.cerebellum_eQTL_fordownload.h5ad 