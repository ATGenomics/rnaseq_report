#!/usr/bin/env bash

SMK_FILE="$HOME/workshop/rnaseq_report/analysis/rnaseq-SE.smk"
if command -v conda &> /dev/null
then
    source activate snakemake
    snakemake -j1 -s ${SMK_FILE} --use-conda
    conda deactivate
else
    echo "Check if conda is installed"
fi


