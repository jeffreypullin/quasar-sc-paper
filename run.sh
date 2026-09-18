#!/bin/bash

export JAVA_HOME=/home/jp2045/software/jdk-23.0.1
export PATH=/home/jp2045/software/jdk-23.0.1/bin:$PATH

mkdir -p data/tmp

./nextflow -log /dev/null run main.nf -resume "$@"
