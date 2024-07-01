# Assess redundancy of the novel CRSs by
# 1. Phylogenetic similarity f contained species
# 2. Overlap in sequence overlaps (relative to genomic coordinates)
# 3. Alignment of the consensus structures with RNAdistance


rule L_redundancy:
    input: 
        'data/K2_motifs.tsv',
        'data/K_motif-tax-pos.tsv',
        'data/J_motif-aln-seqs.tsv.gz',
        'data/J_novel/all_intergenic_regions.tsv.gz',
        # implicit:
        # 'data/F_cmfinder/D_search-seqs/%s/%s'
        # 'data/D_search-seqs/*.fna.gz', filenames but not content
        # 
    output:
        'data/L1_redundancy/no.crs.regions.tsv',
        'data/L1_redundancy/potential.redundant.tsv',
        'data/L1_redundancy/potential.redundant.jaccard.tsv',
        'data/L1_redundancy/relative.overlaps.tsv',
        'data/L1_redundancy/alignment.proportion.tsv',
        'data/L1_redundancy/redundant.candidates.tsv',
        'data/L1_redundancy/redundant.candidates.consensus.tsv',
        'data/L1_redundancy/RNAdistance.pairs.tsv',
        'data/L1_redundancy/RNAdistance.input.txt'
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/L_redundancy.R'
        

# Pairwise RNAdistance on filtered consensus structures
rule L_rnadist:
    input: 
        'data/L1_redundancy/RNAdistance.input.txt'
    output:
        'data/L1_redundancy/RNAdistance.output.txt'
    container: 'viennarna\:2.6.3--py39pl5321h4e691d4_0'
    shell:
      """
      # -B prints alignments for sanitizing
      cat {input} | RNAdistance -B > {output}
      """


# wrap up with results from RNAdistance
        
rule L_redundancy2:
    input: 
        'data/K2_motifs.tsv',
        'data/L1_redundancy/no.crs.regions.tsv',
        'data/L1_redundancy/potential.redundant.tsv',
        'data/L1_redundancy/potential.redundant.jaccard.tsv',
        'data/L1_redundancy/relative.overlaps.tsv',
        'data/L1_redundancy/alignment.proportion.tsv',
        'data/L1_redundancy/redundant.candidates.tsv',
        'data/L1_redundancy/redundant.candidates.consensus.tsv',
        'data/L1_redundancy/RNAdistance.pairs.tsv',
        'data/L1_redundancy/RNAdistance.output.txt'
    output:
        'data/L_rnadistance.tsv',
        'data/L_redundant.tsv',
        'data/L_redundant.jpeg'
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/L_redundancy2.R'
