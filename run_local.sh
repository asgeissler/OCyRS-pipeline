#!/bin/bash
# Helper to run snakemake workflow

# if no parameter passed to script, run all
if [[ $# -eq 0 ]] ; then
    target='all'
else
    target="$@"
fi

# make sure output folder for log files exists
if [ ! -d snakelogs ]; then
  mkdir snakelogs
fi

snakemake                    \
    --use-singularity        \
    --use-conda              \
    --cores all              \
    --resources mem_mb=40000 \
    --keep-going             \
    $target
