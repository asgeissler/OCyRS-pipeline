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


# - Check score distributions
# - Recall of alignment sequences
# - Determine a gathering score
# - List homologs
rule J_homologs:
    input: 
        cmsearch = 'data/J_cmsearch-collected.tsv',
        pos = 'data/J_motif-aln-seq-pos.tsv',
        cats = 'data/J_FDR-categories.tsv'
    output:
        # TODO adapt to new run
        cutoffs = 'data/J_homologs/gathering-scores.tsv',
        homologs = 'data/J_homologs/motif-homologs.tsv',
        fig_jaccard = 'data/J_homologs/overlaps.jpeg',
        fig_boxplot = 'data/J_homologs/score-boxplots.jpeg',
        fig_cor = 'data/J_homologs/score-cor.jpeg',
        fig_powcov = 'data/J_homologs/high-power-covariation.jpeg',
        fig_eval = 'data/J_homologs/eval.jpeg',
        fig_homologs = 'data/J_homologs/motif-homologs.jpeg'
    container: 'renv/renv.sif'
    threads: 16
    conda: 'renv'
    script:
        '../scripts/J_homologs.R'
