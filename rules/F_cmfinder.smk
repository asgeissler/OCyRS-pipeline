# Detect potential motifs in the regions
# 'data/D_search-seqs/{region}.fna.gz'
# and the random background under
# 'data/E_random/{region}.fna.gz'

rule F_cmfinder:
    input:
        'data/{path}/{region}.fna.gz'
    output:
        directory('data/F_cmfinder/{path}/{region}')
    container:
        'cmfinder/cmfinder-0.4.1.18.sif'
    shell:
        """
        PROJ=$(pwd -L)
        # Copy to local tmp and decompress
        tmp=$(mktemp -d)
        cp {input} $tmp/{wildcards.region}.fna.gz
        cd $tmp
        gunzip {wildcards.region}.fna.gz

        # Run CMfinder
        cmfinder04.pl                                                                                                 \
            -skipClustalW -minCandScoreInFinal 10 -combine -fragmentary                                               \
            -commaSepEmFlags x--filter-non-frag,--max-degen-per-hit,2,--max-degen-flanking-nucs,7,--degen-keep,--amaa \
            {wildcards.region}.fna

        # Save to output
        cd $PROJ
        mkdir -p {output}
        # motifs end in numbers
        # (true is needed for strict mode, in case no motifs were detected)
        cp $tmp/*motif*[0-9] {output}/ || true
        # cleanup rest
        rm -rf $tmp

        # set a flag to be absolutely sure that the run was complete
        # (cmfinder is awkward)
        touch {output}/done.flag
        """


# Trigger  CMfinder
rule F_all:
    input:
        lambda wild: D_aggregate(wild, 'data/F_cmfinder/D_search-seqs/{region}')
    output:
        touch('data/F_cmfinder/done.flag')


checkpoint F_collect:
    input:
        'data/F_cmfinder/done.flag'
    output:
        'data/F_cmfinder/motifs.txt'
    shell:
        """
        find data/F_cmfinder/cmfinder -name "*motif*" -exec basename {{}} \; \
            > {output}
        """

