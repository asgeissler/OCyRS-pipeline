# Rules to investigate the space of base-pairing probabilities
# in order to identify potential structures with alternative
# secondary structure conformations


checkpoint M_export:
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
        '../scripts/M_export-msa.R'
    

rule M_alifold:
    input:
        'data/M_alignments-{dir}/{file}.fna'
    output:
        directory('data/M_RNAalifold/{dir}/{file}/')
    container: 'viennarna\:2.6.3--py39pl5321h4e691d4_0'
    shell:
        """
        p=$PWD
        mkdir {output}
        cd {output}
        RNAalifold --input-format=F  --aln --color -p                   \
            --ribosum_scoring                                           \
            --paramFile=/usr/local/share/ViennaRNA/rna_turner2004.par   \
            $p/{input} > stdout.txt
        """


def M_aggregate(wildcards):
    # make sure that the checkpoint finished
    chk = checkpoints.M_export.get().output
    # List files
    xs = glob('data/M_alignments-*/*.fna')
    # Create list of desired RNAalifold output
    return [ 'data/M_RNAalifold/' + \
        i.split('data/M_alignments-')[1][:-len('.fna')] \
        for i in xs ]

rule M_alli:
    input:
        M_aggregate
    output:
        touch('data/M_RNAalifold.done')


