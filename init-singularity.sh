#!/bin/sh
singularity pull https://depot.galaxyproject.org/singularity/busco:5.4.0--pyhdfd78af_2
singularity pull https://depot.galaxyproject.org/singularity/checkm-genome:1.2.0--pyhdfd78af_0
singularity pull https://depot.galaxyproject.org/singularity/cmfinder:0.4.1.9--pl5.22.0_1
singularity pull https://depot.galaxyproject.org/singularity/infernal:1.1.4--pl5321hec16e2b_1
singularity pull https://depot.galaxyproject.org/singularity/muscle:3.8.1551--h7d875b9_6
singularity pull https://depot.galaxyproject.org/singularity/iqtree:2.2.0.3--hb97b32f_0


# Download custom images from ghcr.io 
singularity pull oras://ghcr.io/asgeissler/renv:0.1
singularity pull oras://ghcr.io/asgeissler/treeshrinkenv:0.1
# Restore paths
mv renv_0.1.sif renv/renv.sif
mv treeshrinkenv_0.1.sif treeshrinkenv/treeshrinkenv.sif
