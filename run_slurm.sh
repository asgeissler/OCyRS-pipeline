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

# Get full path to project dir without changing working dir
ABSPATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#echo ${ABSPATH}

snakemake                                                            \
    --use-singularity                                                \
    --jobs 9999 `# limit on jobs in parallel`                        \
    --scheduler greedy `# circumvent stuck on 'select jobs' bug`     \
    --keep-going `# Go on with independent jobs if a job fails.`     \
    --notemp `# required setting by immediate-submit`                \
    --latency-wait 60                                                \
    --immediate-submit `# submit job immediately`                    \
    --cluster-config "config/slurm.yaml"                             \
    `# command to send to batch with dependencies`                   \
    --cluster "sbatch --partition={cluster.partition}                \
        --parsable                                                   \
        \$(if [ -n \"{dependencies}\" ] ; then                       \
            echo \"--dependency=afterok:{dependencies}\" | sed 's/ /,afterok:/g' ;        \
          fi )                                                       \
        --job-name={cluster.job-name}                                \
        --ntasks=1                                                   \
        --cpus-per-task={cluster.cpus}                               \
        --mem={cluster.mem}                                          \
        --time={cluster.time}                                        \
        --output={cluster.output}                                    \
        --error={cluster.error}"                                     \
    $target

