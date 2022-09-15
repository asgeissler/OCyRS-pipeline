# Estimate phylogenetic tree per alignments
rule C_phylo:
    input:
        'data/B_OGs-aln-filtered/{term}.faa.gz'
    output:
        report    = 'data/C_phylo/{term}/report.txt',
        ml        = 'data/C_phylo/{term}/maximumlikelihood.tree',
        mldist    = 'data/C_phylo/{term}/maximumlikelihood.dist',
        splits    = 'data/C_phylo/{term}/bootstrap-splits.nex',
        consensus = 'data/C_phylo/{term}/bootstrap-consensus.tree',
        screen    = 'data/C_phylo/{term}/screen.log'
    container:
        'iqtree\:2.2.0.3--hb97b32f_0'
    resources:
        mem_mb = 20000
    threads: 32
    shell:
        """
        tmp=$(mktemp -d)
        gunzip -c {input} > $tmp/in.faa
        cd $tmp

        iqtree2 -s in.faa -st AA -pre res          \
            -T {threads} -B 1000 -m LG+G -seed 123

        cd {config[project_root]}

        mv $tmp/res.iqtree     {output.report}
        mv $tmp/res.treefile   {output.ml}
        mv $tmp/res.mldist     {output.mldist}
        mv $tmp/res.splits.nex {output.splits}
        mv $tmp/res.contree    {output.consensus}
        mv $tmp/res.log        {output.screen}

        rm -rf $tmp
        """

# Remove outliers from tree
rule C_shrink:
    input:
        'data/C_phylo/{term}/bootstrap-consensus.tree'
    output:
        'data/C_shrink/{term}/output.tree',
        'data/C_shrink/{term}/output_summary.txt',
        'data/C_shrink/{term}/output.txt'
    log: 'snakelogs/C_shrink/{term}.txt'
    container: 'treeshrinkenv/treeshrinkenv.sif'
    conda: 'treeshrinkenv'
    shell:
        """
        run_treeshrink.py                              \
            --tree {input}                             \
            `# Highly recommended for large trees`     \
            --centroid                                 \
            `# After manual inspection, more sensitive`\
            `# filtering than default 0.05 is needed`  \
            --quantiles 0.10                           \
            --outdir 'data/C_shrink/{wildcards.term}/' \
            > {log}
        """

# Compute Pairwise Topology distances and check out MDS plot
# incl. comparison with reference trees.
# Challenge: Passing tree files as parameters will result in a too large file
# for sbatch to handle.
# => Implicit dependency via flag files

rule C_raw_flag:
    input:
        lambda wild: B_aggregate_OG(wild, 'data/C_phylo/{term}/bootstrap-consensus.tree')
    output:
        touch('data/C_phylo/flag.done')


rule C_shrink_flag:
    input:
        lambda wild: B_aggregate_OG(wild, 'data/C_shrink/{term}/output.txt'),
    output:
        touch('data/C_shrink/flag.done')


rule C_space:
    input:
        'data/C_phylo/flag.done',
        'data/C_shrink/flag.done'
    output:
        raw = 'data/C_space/pairwise-distances-raw.tsv',
        shrunk = 'data/C_space/pairwise-distances-shrunk.tsv',
        mds = 'data/C_space/mds-data.tsv',
        mdsfig = 'data/C_space/mds.jpeg'
    log: 'snakelogs/C_space.txt'
    container: 'renv/renv.sif'
    threads: 32
    conda: 'renv'
    script:
        '../scripts/C_space.R'

