# OCyRS: Orthologues cyanobacterial gene-associated RNA structure motifs

This repository corresponds to the pipline of the scientific publication:


A cyanobacteria phylum-wide genomic screen for RNA structure motifs
adjacent to orthologues genes  
AS  Geissler, EC Alvarez, NU Frigaard, J Gorodkin, and SE Seemann  
*in preparation*


## Purpose

This pipline download public cyanobacterial genome sequence and annotations
from
[proGenomes](https://progenomes.embl.de/).
Anchored to orthologues genes, potential adjacently located
RNA structure motif are searched for with
[CMfinder](http://bio.cs.washington.edu/yzizhen/CMfinder/).
The pipline includes additional steps to score the motifs and check
for their phylogenetic/biological plausibility with
[R-scape](http://eddylab.org/R-scape/).
As such, the motifs are compared with known bacterial
[RNA families](https://rfam.xfam.org/),
the 
[Infernal](http://eddylab.org/infernal/)
predicted genomic location of these families,
and the respective scores.


## Requirements

This pipeline is implemented in a reproducible
[snakemake](https://snakemake.github.io/)
workflow, which automatically handles the further
software dependencies. The only hard-requirement are

1. A powerful machine or cluster (at least 20 CPU cores 20+ GB RAM)
2. Singularity version $\ge 3.9$
3. [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) + Mamba

### MacOS

To install singularity, please consider using
a [Vagrant](https://www.vagrantup.com/) box for
virtualisation.
The installation is straight-forward with 
[Homebrew](https://brew.sh/):

        $ brew install --cask virtualbox vagrant vagrant-manager

With the information in `Vagrantfile`, singularity can be installed and
activated in the current working directory with:

        $ vagrant up
        $ vagrant ssh
        $ cd /vagrant/

Finally, the conda environment must be installed within the virtual box.

        $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
        $ sh Miniconda3-latest-Linux-x86_64.sh

Due to the requirement 'For changes to take effect, close and re-open your current shell.', you must exit and reopen the shell before being able to install mamba.

        $ exit
        $ vagrant ssh
        $ cd /vagrant/
        $ conda install -n base -c conda-forge mamba

Afterward, the virtual box is fully setup to work with the pipeline.

## Dependencies

Download all further software dependencies with:

        $ sh init-singularity.sh

The script downloads the following software images from the
[Galaxy project](https://galaxyproject.org/)
software repository:

1. BUSCO v$5.4.0$
1. CMfinder v$0.4.1.9$
1. Infernal v$1.1.4$
1. Snakemake v$7.9.0$





