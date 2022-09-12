This folder described the environment of the R packages. To ensure maximal reproducibility, the 'raw' description file was locked in place on 29th July 2022 with

        conda lock -p linux-64 -f renv.yml

The resulting file `conda-lock.yml` can then be used to build a singularity
image, as described in `Singularity`.

        sudo singularity build renv.sif Singularity

A pre-build version of this environment is attached to the github 
packages [repository](https://github.com/asgeissler/OCyRS-pipeline).

