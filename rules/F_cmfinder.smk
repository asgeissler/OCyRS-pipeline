# Detect potential motifs in the regions
# 'data/D_search-seqs/{region}.fna.gz'
# and the random background under
# 'data/E_search-shuffled/{region}.fna.gz'

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
            {wildcards.region}.fna     || true

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


# Trigger  CMfinder for biological sequneces
rule F_bioruns:
    input:
        lambda wild: D_aggregate(wild, 'data/F_cmfinder/D_search-seqs/{region}')
    output:
        'data/F_cmfinder/D_search-seqs/motifs.txt'
    shell:
        """
        find data/F_cmfinder/D_search-seqs -name "*motif*" > {output}
        """

# Trigger CMfinder for the shuffeled sequences
rule F_random:
    input:
        lambda wild: E_bg_models(
            wild,
            'data/F_cmfinder/E_search-shuffled_{seed}/{{region}}'.format(**wild)
        )
    output:
        'data/F_cmfinder/E_search-shuffled_{seed}/motifs.txt'
    shell:
        """
        find data/F_cmfinder/E_search-shuffled_{wildcards.seed}  \
           -name "*motif*"  > {output}
        """

checkpoint F_demerge:
    input:
        'data/F_cmfinder/D_search-seqs/motifs.txt',
        *['data/F_cmfinder/E_search-shuffled_{}/motifs.txt'.format(i) \
          for i in E_seeds]
    output:
        'data/F_cmfinder/motifs-demerged.tsv'
    log: 'snakelogs/F_demerged.txt'
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/F_demerge.R'


# allow access to fields in demerged.tsv
# - path
# - filename
# - dir (eg D_search-seqs)
# - region
def F_aggregate(wildcards, x):
    chk = checkpoints.F_demerge.get().output
    df = pd.read_csv('data/F_cmfinder/motifs-demerged.tsv', sep = '\t')
    return sample_wise(x, df)

