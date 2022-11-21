rule J_collect:
    input: 
        'data/I_cmsearch/runs.done'
    output:
        'data/J_cmsearch-collected.tsv',
    container: 'renv/renv.sif'
    threads: 32
    conda: 'renv'
    script:
        '../scripts/J_collect.R'



rule J_calibrate:
    input: 
        collected = 'data/J_cmsearch-collected.tsv',
        motifs = 'data/F_cmfinder/motifs-demerged.tsv'
    output:
        motifpos = 'data/J_motif-aln-seq-pos.tsv',
        cutoffs = 'data/J_score-cutoffs.tsv',
        homologs = 'data/J_motif-homologs.tsv',
        fig = 'data/J_scoring.png'
    container: 'renv/renv.sif'
    threads: 16
    conda: 'renv'
    script:
        '../scripts/J_check-aln.R'
