#!/bin/bash
# Helper to run snakemake workflow

# if no parameter passed to script, run all
if [[ $# -eq 0 ]] ; then
    target='all'
else
    target="$@"
fi

singularity exec                  \
    --bind /vagrant:/vagrant      \
    --pwd /vagrant                \
    snakemake:7.9.0--hdfd78af_0   \
    snakemake                     \
    --use-singularity             \
    --cores all                   \
    $target
