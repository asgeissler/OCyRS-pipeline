# Estimate tree per alignments
rule C_phylo:
    input:
        'data/B_OGs-aln/{term}.faa.gz'
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
    threads: 16
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

# trigger tree per alignment, using helper function from B_*
rule C_aggregate:
    input:
        lambda wild: B_aggregate_OG(wild, 'data/C_phylo/{term}/report.txt')
    output:
        touch('data/C_aggregate.flag')

