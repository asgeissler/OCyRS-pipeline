#!/bin/bash
# Helper to run snakemake workflow

# if no parameter passed to script, run all
if [[ $# -eq 0 ]] ; then
    target='all'
else
    target="$@"
fi

# Get full path to project dir without changing working dir
ABSPATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#echo ${ABSPATH}

snakemake                                                 \
    --use-singularity                                     \
    `# ensure containers working dir binds to project`    \
    --singularity-args "--bind $ABSPATH:/mnt --pwd /mnt"  \
    --cores all                                           \
    $target
