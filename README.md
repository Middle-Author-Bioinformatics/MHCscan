# MHCscan
Software tool for identification of multiheme cytochromes

### [Install Miniconda](https://docs.conda.io/en/latest/miniconda.html)
### Installation
    git clone https://github.com/Arkadiy-Garber/MHCscan/tree/main
    conda create -n mhcscan -c bioconda -c conda-forge -c defaults hhsuite prodigal --yes
    
#### *Put MHCscan main repo and bin into your PATH variable*

    export PATH="/path/to/MHCscan:$PATH"
    export PATH="/path/to/MHCscan/bin:$PATH"

*(recommended to add this to your .bashrc, .bash_profile, .profile, or .zshrc scripts)*

### [Download databases for hhsuite](https://github.com/soedinglab/hh-suite)
We recommend using PDB70, which can be downloaded with the following command (make sure you have > 60 Gb of free storage space)

    wget https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_latest.tar.gz

### Usage (with annotated genome)
    MHCscan.sh -i proteins.faa -g annotation.gff -o output.csv -t 16 --prot

### Usage (with unannotated genome)
    MHCscan.sh -i contigs.fasta -o output.csv -t 16

![pipeline](https://github.com/Arkadiy-Garber/MHCscan/blob/main/pipeline.png)


## Gram-positive candidate EET wrapper

A new `mhc-eet` CLI combines MHCscan sequence characterization, FeGenie family
HMMs, FEET motif compatibility, and coordinate-based surface-domain neighborhoods.
It keeps same-protein domains, nearby anchoring machinery and optional localization
evidence separate. The default broad screen is at least two canonical CXXCH motifs.

See [installation, database setup, commands and interpretation](docs/GRAM_POSITIVE_EET.md),
[repository/license audit](docs/REPOSITORY_AUDIT.md), and [validation](docs/VALIDATION.md).
This extension was added in November 2026; the original command remains available.
