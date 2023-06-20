# MHCscan
Software tool for identification of multiheme cytochromes

### [Install Miniconda](https://docs.conda.io/en/latest/miniconda.html)
### Installation
    conda create -n mhcscan -c bioconda -c conda-forge -c defaults hhsuite prodigal --yes
### [Download databases for hhsuite](https://github.com/soedinglab/hh-suite)

### Usage
    MHCscan.sh -i proteins.faa -g annotation.gff -o output.csv -t 16 --prot

![pipeline](https://github.com/Arkadiy-Garber/MHCscan/blob/main/pipeline.png)

