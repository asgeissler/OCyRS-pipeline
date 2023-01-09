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


# - Check score distributions
# - Recall of alignment sequences
# - Determine a gathering score
# - List homologs
rule J_homologs:
    input: 
        cmsearch = 'data/J_cmsearch-collected.tsv',
        pos = 'data/J_motif-aln-seq-pos.tsv',
        fdr = 'data/I_fdr.tsv',
        scores = 'data/H_scores.tsv',
        cmstat = 'data/I_cmstat.tsv'
    output:
        cutoffs = 'data/J_gathering-scores.tsv',
        homologs = 'data/J_motif-homologs.tsv',
        fig_jaccard = 'data/J_overlaps.jpeg',
        fig_boxplot = 'data/J_score-boxplots.jpeg',
        fig_cor = 'data/J_score-cor.jpeg',
        fig_powcov = 'data/J_high-power-covariation.jpeg',
        fig_eval = 'data/J_eval.jpeg',
        fig_homologs = 'data/J_motif-homologs.jpeg'
    container: 'renv/renv.sif'
    threads: 16
    conda: 'renv'
    script:
        '../scripts/J_homologs.R'
