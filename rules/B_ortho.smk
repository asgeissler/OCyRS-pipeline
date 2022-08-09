
rule B_cand:
    input:
        tax = 'data/A_representatives/taxonomy.tsv',
        genes = 'data/A_representatives/genes.tsv.gz',
        kegg = 'data/A_representatives/kegg.tsv.gz'
    output:
        fig = 'data/B_OGs.jpeg',
        tbl = 'data/B_OGs.tsv'
    log: 'snakelogs/B_cand.txt'
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/B_orthology-groups.R'

