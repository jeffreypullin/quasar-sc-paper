#!/usr/bin/env bash

find /home/jp2045/quasar-sc-paper/data/nextflow -mindepth 1 -maxdepth 1 ! -name conda -exec rm -rf {} +
