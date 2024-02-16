#!/usr/bin/env bash

# install miniconda in bash:

# This script will install miniconda and it will
# create a snakemake and nextflow environments.

# the custom installation use 'source' instead conda to activate environments

export PATH=$HOME/bin/miniconda3/bin:$PATH

if command -v conda &> /dev/null
then
    conda install -n base -c conda-forge mamba
    conda create -n snakemake -c conda-forge -c bioconda snakemake
else
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/bin/miniconda3
    rm Miniconda3-latest-Linux-x86_64.sh

cat <<"EOF" >> $HOME/.bashrc
export PATH=$HOME/bin/miniconda3/bin:$PATH
EOF

fi