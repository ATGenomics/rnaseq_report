#!/usr/bin/env bash

source activate snakemake

snakemake -j1 -s rnaseq-SE.smk --use-conda

conda deactivate