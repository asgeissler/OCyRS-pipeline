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


#rule J_calibrate:
#    input: 
#        collected = 'data/J_cmsearch-collected.tsv',
#        pos = 'data/J_motif-aln-seq-pos.tsv'
#    output:
#        cutoffs = 'data/J_score-cutoffs.tsv',
#        homologs = 'data/J_motif-homologs.tsv',
#        fig = 'data/J_scoring.png'
#    container: 'renv/renv.sif'
#    threads: 16
#    conda: 'renv'
#    script:
#        '../scripts/J_check-aln.R'
