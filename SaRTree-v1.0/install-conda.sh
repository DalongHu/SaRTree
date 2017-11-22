#!/bin/bash
echo "SaRTree installation strating...\n"
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda config --add channels esteinig
conda create -y --name env_sartree perl perl-yaml-libyaml perl-bioperl perl-findbin perl-getopt-long raxml mauve beast
echo "installation finished, please use 'source activate env_sartree' to access SaRTree environment\n"
