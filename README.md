# MHCscan
Software tool for identification of multiheme cytochromes

### [Install Miniconda](https://docs.conda.io/en/latest/miniconda.html)
### Installation
    git clone https://github.com/Arkadiy-Garber/MHCscan/tree/main
    conda create -n mhcscan -c bioconda -c conda-forge -c defaults hhsuite prodigal --yes
    
#### *Put YfGenie main repo and bin into your PATH variable (recommended to add this to your .bashrc, .bash_profile, .profile, or .zshrc scripts)*
    export PATH="/path/to/MHCscan:$PATH"
    export PATH="/path/to/MHCscan/bin:$PATH"
    
### [Download databases for hhsuite](https://github.com/soedinglab/hh-suite)
We recommend using PDB70, which can be downloaded with the following command (make sure you have > 60 Gb of free storage space)

    wget https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_latest.tar.gz

### Usage
    MHCscan.sh -i proteins.faa -g annotation.gff -o output.csv -t 16 --prot

![pipeline](https://github.com/Arkadiy-Garber/MHCscan/blob/main/pipeline.png)

