
# Determine the candidate OGs
rule B_cand:
    input:
        tax = 'data/A_representatives/taxonomy.tsv',
        genes = 'data/A_representatives/genes.tsv.gz',
        kegg = 'data/A_representatives/kegg.tsv.gz'
    output:
        fig = 'data/B_OGs.jpeg',
        tbl = 'data/B_OGs.tsv'
    log: 'snakelogs/B_cand.txt'
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/B_orthology-groups.R'

# Write out sequences of candidate OGs
checkpoint B_cand_seqs:
    input:
        tbl = 'data/B_OGs.tsv',
        seqs = 'data/A_representatives'
    output:
        directory('data/B_OGs')
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/B_write-seqs.R'


# align faa sequences
rule B_aln:
    input:
        'data/B_OGs/{term}.faa.gz'
    output:
        'data/B_OGs-aln/{term}.faa.gz'
    container:
        'muscle\:3.8.1551--h7d875b9_6'
    shell:
        """
        tmp=$(mktemp)
        gunzip -c {input} > $tmp
        muscle -in $tmp -out $tmp.out
        gzip $tmp.out
        mv $tmp.out.gz {output}
        """

def B_aggregate_OG(wildcards, x):

    """
    List for all OGs the fasta sequence file with
        x = 'data/B_OGs/{term}.faa.gz'):
    """
    # make sure that the checkpoint finished
    chk = checkpoints.B_cand_seqs.get().output[0]
    # load taxonomy table without dots in names
    df = pd.read_csv('data/B_OGs.tsv', sep = '\t')
    df.columns = df.columns.str.replace('.', '_')
    # list genome files
    return sample_wise(x, df)



rule B_aln_done:
    input:
        lambda wild: B_aggregate_OG(wild, 'data/B_OGs-aln/{term}.faa.gz')
    output:
        touch('data/B_OGs-aln/done.flag')
