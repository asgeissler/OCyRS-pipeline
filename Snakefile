import pandas as pd
import os
from glob import glob

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
include: 'rules/D_search-seqs.smk'
include: 'rules/E_background.smk'
include: 'rules/F_cmfinder.smk'
include: 'rules/G_rfam.smk'
include: 'rules/G2_rnie.smk'
include: 'rules/H_scores.smk'
include: 'rules/I_fdr.smk'
include: 'rules/J_recall.smk'

# Target rules
# Initial proGenomes download is separate to keep remaining rules streamlined
rule download:
    input:
        'data/A_progenomes/genomes.fna.gz',
        'data/A_progenomes/representatives.txt',
        'data/A_representatives/'

rule all:
    input:
        # initiate extra genome quality check
        'data/A_checkm/checkm_summary.tsv',
        # extract OG, their sequences, and align
        'data/B_OGs.tsv',
        'data/B_OGs-aln-filtered/done.flag',
        'data/B_seqids/seqid.tsv',
        # build trees and compare their topology distances
        'data/C_space/pcoa.jpeg',
        # export adjacent sequences
        'data/D_intergenic.jpeg',
        'data/D_search-seqs-aln/done.flag',
        # make random background
        'data/E_search-shuffled.done-flag',
        'data/E_search-shuffled_stat.flag',
        # CMfinder runs
        'data/F_cmfinder/motifs-demerged.tsv',
        # Also make a bacterial Rfam screen in all genomes
        'data/G_rfam-bacteria/download.done',
        'data/G_rfam-bacteria-seeds/download.done',
        'data/G_rfam-cmsearch/runs.done',
        'data/G_rfam-cmsearch.tsv.gz',
        # Prediction of terminators
        'data/G2_terminators/done.flag',
        # Score sto files
        'data/H_scores/done.flag',
        'data/H_scores.tsv',
        #  Estimate FDRs, and scan genomes for candidates
        'data/I_fdr.tsv',
        'data/I_candidate-models/done.cm',
        'data/I_cmsearch/runs.done',
        # Calibrate score cutoff per model
        # then recall rfam hits
        'data/J_cmsearch-collected.tsv',
        'data/J_scoring.png'



