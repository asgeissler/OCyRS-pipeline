# OCyRS: Orthologues cyanobacterial gene-associated RNA structure motifs

This repository corresponds to the pipeline of the scientific publication:


A cyanobacteria phylum-wide genomic screen for RNA structure motifs
adjacent to orthologues genes  
AS  Geissler, EC Alvarez, NU Frigaard, J Gorodkin, and SE Seemann  
*in preparation*


## Purpose

This pipeline download public cyanobacterial genome sequence and annotations
from
[proGenomes](https://progenomes.embl.de/).
Anchored to orthologues genes, potential adjacently located
RNA structure motif are searched for with
[CMfinder](http://bio.cs.washington.edu/yzizhen/CMfinder/).
The pipeline includes additional steps to score the motifs and check
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
3. [Snakemake](https://snakemake.github.io/) $\ge 7.9.0$
4. Git [large file support](https://git-lfs.github.com/) plugin

### Linux

Please use either the binary packages coming with your distribution's package
system or follow the installation procedures listed by the respective
project websites.

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
        $ mamba install -n base -c bioconda -c conda-forge snakemake

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


A container with various R packages `renv/renv.sif`
is pre-build as part of this pipeline.
The corresponding conda environment was for maximal reproducibility
[conda-lock'ed](https://github.com/conda-incubator/conda-lock),
such that the exact versions are specified for linux `x64`
computer architectures, in case you would like to re-build the container.



## How to use 

***Before you start:***
Please adjust the `project_root` variable in the `config/config.yaml` file to
your project path.


There are currently two helper script to aid starting snakemake. Use either

1. `$ bash run_local.sh` for running snakemake on a single powerful machine.
2. `$ bash run_slurm.sh` to submit jobs to a slurm cluster.

In case of the latter, also consider adjusting `config/slurm.yaml` to your
environment (*e.g.* partition name). To adapt a similar helper script for
alternative cluster solutions, please refer to the
[snakemake handbook](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).


