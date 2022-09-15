
checkpoint D_seqs:
    input:
        'data/A_representatives/genes.tsv.gz',
        # depend on flag to prevent too long dependency list for sbatch
        'data/C_shrink/flag.done'
    output:
        seqs = directory('data/D_search-seqs/'),
        intergenic = 'data/D_intergenic.jpeg'
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/D_search-seqs.R'

# Similar to  `B_aggregate_OG` wait for search sequence generation checkpoint
# to fill sequence
# make sure alignment runs after the B_cand_seqs checkpoint
def D_aggregate(wildcards, x):
    # make sure that the checkpoint finished
    chk = checkpoints.D_seqs.get().output
    # Check which regions where exported form filesystem
    xs = glob('data/D_search-seqs/*fna.gz')
    xs = [ os.path.basename(i) for i in xs ]
    xs = [ i.split('.')[0] for i in xs ]
    # expand desired string
    df = pd.DataFrame({'region': xs})
    return sample_wise(x, df)


# similar to B_align rule
rule D_align:
    input:
        'data/D_search-seqs/{region}.fna.gz'
    output:
        'data/D_search-seqs-aln/{region}.fna.gz'
    container:
        'muscle\:3.8.1551--h7d875b9_6'
    shell:
        """
        tmp=$(mktemp)
        gunzip -c {input} > $tmp
        muscle -in $tmp -out $tmp.out -maxiters 2
        gzip $tmp.out
        mv $tmp.out.gz {output}
        """


rule D_collect:
    input:
        lambda wild: D_aggregate(wild, 'data/D_search-seqs-aln/{region}.fna.gz')
    output:
        touch('data/D_search-seqs-aln/done.flag')

