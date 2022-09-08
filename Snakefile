import pandas as pd
import os

from snakemake.utils import validate, min_version
min_version("7.9.0")

# Load configuration
configfile: 'config/config.yaml'
validate(config, 'schemas/config.schema.yaml')


# Helper to write output/input of rules
def sample_wise(x, df):
    "Expand wildcards in x as a generator for data in df"
    return list(df.apply(
        lambda row: x.format(**row),
        axis = 1

    ))

include: 'rules/A_genomes.smk'
include: 'rules/B_ortho.smk'
include: 'rules/C_phylo.smk'

# Target rules
rule all:
    input:
        'data/A_progenomes/genomes.fna.gz',
        'data/A_progenomes/representatives.txt',
        'data/A_representatives/',
        #
        ##'data/A_checkm/checkm_summary.tsv',
        #
        'data/B_OGs.tsv',
        'data/B_OGs-aln-filtered/done.flag',
        'data/B_seqids/seqid.tsv',
        'data/C_aggregate.flag'
