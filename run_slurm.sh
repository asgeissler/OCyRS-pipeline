#!/bin/bash
# Helper to run snakemake workflow

# if no parameter passed to script, run all
if [[ $# -eq 0 ]] ; then
    target='all'
else
    target="$@"
fi

echo "Targets: $target"

# make sure output folder for log files exists
if [ ! -d slurmlogs ]; then
  mkdir -p slurmlogs
fi

snakemake                                                            \
    --use-singularity                                                \
    --jobs 999 `# limit on jobs in parallel`                         \
    --scheduler greedy `# circumvent stuck on 'select jobs' bug`     \
    --keep-going `# Go on with independent jobs if a job fails.`     \
    --latency-wait 60                                                \
    --cluster-config "config/slurm.yaml"                             \
    `# command to send to batch with dependencies`                   \
    --cluster "sbatch --partition={cluster.partition}                \
        --job-name={cluster.job-name}                                \
        --ntasks=1                                                   \
        --cpus-per-task={cluster.cpus}                               \
        --mem={cluster.mem}                                          \
        --time={cluster.time}                                        \
        --output={cluster.output}                                    \
        --error={cluster.error}"                                     \
    $target

