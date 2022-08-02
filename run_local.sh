#!/bin/bash
# Helper to run snakemake workflow

# if no parameter passed to script, run all
if [[ $# -eq 0 ]] ; then
    target='all'
else
    target="$@"
fi

snakemake                                                 \
    --use-singularity                                     \
    --cores all                                           \
    $target
