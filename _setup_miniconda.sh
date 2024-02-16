#!/usr/bin/env bash

# install miniconda in bash:

# This script will install miniconda and it will
# create a snakemake and nextflow environments.

# the custom installation use 'source' instead conda to activate environments
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/bin/miniconda3
rm Miniconda3-latest-Linux-x86_64.sh

cat <<"EOF" >> $HOME/.bashrc
export PATH=$HOME/bin/miniconda3/bin:$PATH
EOF

export PATH=$HOME/bin/miniconda3/bin:$PATH

# conda create -n nextflow -yc bioconda nextflow

conda create -n snakemake -c conda-forge -c bioconda snakemake