#!/usr/bin/env bash

printf "\n    ${GREEN}Setting up conda environment...${NC}\n\n"

## adding conda channels
conda config --add channels defaults 2> /dev/null
conda config --add channels bioconda 2> /dev/null
conda config --add channels conda-forge 2> /dev/null
conda config --add channels au-eoed 2> /dev/null

conda create -n mhcscan -c bioconda -c conda-forge -c defaults hhsuite prodigal --yes

## activating environment
conda activate mhcscan

## creating directory for conda-env-specific source files
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d

## adding paths:
echo '#!/bin/sh'" \
export PATH=\"$(pwd):"'$PATH'\"" \
export PATH=\"$(pwd)/bin\"" >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

# re-activating environment so variable and PATH changes take effect
conda activate mhcscan


printf "\n        ${GREEN}DONE!${NC}\n\n"
