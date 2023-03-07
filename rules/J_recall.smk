rule J_categorize:
    input:
        fdr = 'data/I_fdr.tsv',
        scores = 'data/H_scores.tsv',
        cmstat = 'data/I_cmstat.tsv',
        rfam = 'data/G_rfam-cmsearch.tsv.gz'
    output:
        cor = 'data/J_FDR-score-correlations.pdf',
        cut = 'data/J_ECDF-covariation-power.png',
        dist = 'data/J_category-scores.png',
        cats = 'data/J_FDR-categories.tsv'
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/J_categorize.R'
    

rule J_seq_pos:
    input: 
        'data/I_fdr.tsv'
    output:
        pos = 'data/J_motif-aln-seq-pos.tsv',
        seq = 'data/J_motif-aln-seqs.tsv.gz'
    container: 'renv/renv.sif'
    threads: 16
    conda: 'renv'
    script:
        '../scripts/J_aln-seq-pos.R'


