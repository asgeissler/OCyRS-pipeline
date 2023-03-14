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


rule J_novel:
    input: 
        'data/G_rfam-cmsearch.tsv.gz',
        'data/G2_terminators.tsv.gz',
        'data/J_motif-aln-seq-pos.tsv',
        'data/J_FDR-categories.tsv',
        # additionally links to the paths of the
        # genes, genomes, and trees
        # are implicit, since there should be no way to get to the FDRs without
        # them
        # 'data/A_representatives/genes.tsv.gz'
        # 'data/A_representatives/*/genome.fna.gz'
        # 'data/C_shrink/*/output.tree'
    output:
        'data/J_novel/all_intergenic_regions.tsv.gz',
        'data/J_novel/references_inside-of_intergenic_regions.tsv.gz',
        'data/J_novel/reference-motif-overlaps.jpg',
        'data/J_novel/reference-motif-overlap-stats.tsv',
        'data/J_novel/recall-precision-plot.jpg',
        'data/J_novel/overview.tsv',
        'data/J_novel/potentially-novel-motifs.tsv'
    container: 'renv/renv.sif'
    threads: 16
    conda: 'renv'
    script:
        '../scripts/J_recall.R'

