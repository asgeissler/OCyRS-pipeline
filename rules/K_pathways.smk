# Associate potentially novel motifs to pathways as species information

rule K_pathways:
    input: 
        'data/J_novel/potentially-novel-motifs.tsv',
        'data/J_novel/references_inside-of_intergenic_regions.tsv.gz',
        'data/J_motif-aln-seq-pos.tsv',
        'data/J_FDR-categories.tsv',
        'data/A_representatives/taxonomy.tsv'
    output:
        'data/K_ko-path.tsv',
        'data/K_motif-path.tsv',
        'data/K_motif-tax-pos.tsv',
        'data/K_motif-tax.tsv',
        'data/K_overview.tsv'
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/K_pathways.R'