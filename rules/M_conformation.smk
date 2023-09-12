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
        aln = 'data/M_alignments-{dir}/{motif}/aln.fna',
        struct = 'data/M_alignments-{dir}/{motif}/structure.txt'
        # optional:
        #'data/M_alignments-{dir}/{motif}/tree.txt'
    output:
        rel = 'data/M_PETfold/{dir}/{motif}/reliabilities.txt',
        std = 'data/M_PETfold/{dir}/{motif}/output.txt'
    container: 'petfold-2.2.simg'
    shell:
        """
        X='data/M_alignments-{wildcards.dir}/{wildcards.motif}/tree.txt'
        # Optional check for Tree (not available for Rfam)
        if [[ -f "$X" ]] ; then
            opt="--settree $X"
        else
            opt=''
        fi

        PETfold $opt                  \
            --setstruc {input.struct} \
            --fasta {input.aln}       \
            --ppfile {output.rel}     \
            > {output.std}

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


