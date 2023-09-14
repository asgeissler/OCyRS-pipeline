#!/bin/sh
singularity pull https://depot.galaxyproject.org/singularity/checkm-genome:1.2.0--pyhdfd78af_0
singularity pull https://depot.galaxyproject.org/singularity/infernal:1.1.4--pl5321hec16e2b_1
singularity pull https://depot.galaxyproject.org/singularity/muscle:3.8.1551--h7d875b9_6
singularity pull https://depot.galaxyproject.org/singularity/iqtree:2.2.0.3--hb97b32f_0
singularity pull https://depot.galaxyproject.org/singularity/viennarna:2.6.3--py39pl5321h4e691d4_0


# Download custom images from ghcr.io 
singularity pull oras://ghcr.io/asgeissler/renv:0.2
singularity pull oras://ghcr.io/asgeissler/treeshrinkenv:0.1
singularity pull oras://ghcr.io/asgeissler/cmfinder:0.4.1.18
singularity pull oras://ghcr.io/asgeissler/sissiz:3.0
singularity pull oras://ghcr.io/asgeissler/rscape:2.0.0.j
singularity pull oras://ghcr.io/asgeissler/petfold:2.2

# Restore paths to where the rules expect the images to be
mv renv_0.2.sif renv/renv.sif
mv treeshrinkenv_0.1.sif treeshrinkenv/treeshrinkenv.sif
mv cmfinder_0.4.1.18.sif cmfinder/cmfinder-0.4.1.18.sif
mv sissiz_3.0.sif sissiz/sissiz-3.0.sif
mv rscape_2.0.0.j.sif rscape/rscape-2.0.0.j.sif
mv petfold_2.2.sif petfold-2.2.simg
