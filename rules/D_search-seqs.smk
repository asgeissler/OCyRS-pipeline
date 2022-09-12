
rule D_seqs:
    input:
        'data/A_representatives/genes.tsv.gz',
        # depend on flag to prevent too long dependency list for sbatch
        'data/C_shrink/flag.done'
    output:
        seqs = lambda wild: B_aggregate_OG(wild, 'data/C_search-seqs/{term}.fna.gz'),
        intergenic = 'data/D_intergenic.jpeg'
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/D_search-seqs.R'
