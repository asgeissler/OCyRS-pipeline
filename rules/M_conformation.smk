# Rules to investigate the space of base-pairing probabilities
# in order to identify potential structures with alternative
# secondary structure conformations


checkpoint M_prepare:
    input:
        'data/J_novel/potentially-novel-motifs.tsv',
        'data/G_rfam-cmsearch.tsv.gz'
    output:
        directory('data/M_alignments-candidates/'),
        directory('data/M_alignments-Rfam/')
    container: 'renv/renv.sif'
    threads: 16
    conda: 'renv'
    script:
        '../scripts/M_prepare.R'
    


rule M_petfold:
    input:
        'data/M_alignments-{dir}/{motif}/aln.fna'
    output:
        rel = 'data/M_PETfold/{dir}/{motif}/reliabilities.txt',
        std = 'data/M_PETfold/{dir}/{motif}/output.txt'
    container: 'petfold-2.2.simg'
    shell:
        """
        (PETfold --fasta {input} --ppfile {output.rel} > {output.std}) || true
        # graceful fail for the Rfam families with < 3 entires
        touch {output.rel}
        touch {output.std}
        """


def M_aggregate(wildcards):
    # make sure that the checkpoint finished
    chk = checkpoints.M_prepare.get().output
    # List files
    xs = glob('data/M_alignments-*/*')
    # Create list of desired PETfold output
    return [ 'data/M_PETfold/' +                               \
        os.path.dirname(i).split('data/M_alignments-')[1] +    \
        '/' +                                                  \
        os.path.basename(i) +                                  \
        '/output.txt'                                          \
        for i in xs               ]

rule M_alli:
    input:
        M_aggregate
    output:
        touch('data/M_PETfold.done')


rule M_collect:
    input:
        'data/M_PETfold.done'
    output:
        'data/M_PETfold-stats.png',
        'data/M_PETfold.tsv'
    container: 'renv/renv.sif'
    threads: 16
    conda: 'renv'
    script:
        '../scripts/M_petfold-overview.R'

