
rule I_fdr:
    input:
        'data/H_scores/done.flag',
        'data/H_scores.tsv'
    output:
        overview = 'data/I_overview.png',
        est = 'data/I_fdr.png',
        fdr = 'data/I_fdr.tsv',
        stat = 'data/I_alignment-stats.png'
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/I_fdr.R'

