
rule A_download:
    output:
        protected('data/A_progenomes/genomes.fna.gz'),
        protected('data/A_progenomes/genes.fna.gz'),
        protected('data/A_progenomes/proteins.faa.gz'),
        protected('data/A_progenomes/genes.tsv'),
        protected('data/A_progenomes/groups.tsv'),
        protected('data/A_progenomes/antibiotic.tsv'),
        protected('data/A_progenomes/specI_lineageNCBI.tab'),
        protected('data/A_progenomes/specI_clustering.tab'),
        protected('data/A_progenomes/proteins.representatives.fasta.gz')
    shell:
        """
        mkdir -p data/A_progenomes
        cd data/A_progenomes

        # lookup table of download mode to file name
        typeset -A names
        names=(
            [ca]='genomes.fna.gz'
            [ga]='genes.fna.gz'
            [pa]='proteins.faa.gz'
        )
        for k in "${{!names[@]}}"; do
            wget \
                "https://progenomes.embl.de/dumpSequence.cgi?p={config[phylum]}&t=$k&a=phylum" \
                -O ${{names[$k]}}
        done

        # repeat for annotations
        typeset -A names2
        names2=(
            [ga]='genes.tsv'
            [ae]='groups.tsv'
            [aa]='antibiotic.tsv'
        )
        for k in "${{!names2[@]}}"; do
            wget \
                "https://progenomes.embl.de/dumpAnnotation.cgi?p={config[phylum]}&t=$k&a=phylum" \
                -O ${{names2[$k]}}
        done

        # download remaining meta info
        wget 'https://progenomes.embl.de/data/proGenomes2.1_specI_lineageNCBI.tab' -O specI_lineageNCBI.tab
        wget 'https://progenomes.embl.de/data/proGenomes2.1_specI_clustering.tab' -O specI_clustering.tab

        # needed to determine what the represantitive genomes are
        wget 'https://progenomes.embl.de/data/repGenomes/freeze12.proteins.representatives.fasta.gz' -O proteins.representatives.fasta.gz
        """

rule A_replist:
    input:
        'data/A_progenomes/proteins.representatives.fasta.gz'
    output:
        protected('data/A_progenomes/representatives.txt')
    shell:
        """
        # extract the taxid.bioproject identifiers for representative genomes
        foo=$(mktemp)
        gunzip -c {input} > $foo
        grep '^>' $foo > $foo.deflines
        sed "s,^>,,g" $foo.deflines | \
            sed "s,\t.*$,,g"        | \
            sed "s,\.[^.]*$,,g"     | \
            uniq  > $foo.near

        sort $foo.near | uniq > $foo.done
        cp $foo.done {output}

        #rm $foo*
        """

checkpoint A_representatives:
    input:
        rep = 'data/A_progenomes/representatives.txt',
        tax = 'data/A_progenomes/specI_lineageNCBI.tab',
        genomes_seq = 'data/A_progenomes/genomes.fna.gz',
        proteins_seq = 'data/A_progenomes/proteins.faa.gz',
        genes_seq = 'data/A_progenomes/genes.fna.gz',
        genes = 'data/A_progenomes/genes.tsv',
        groups = 'data/A_progenomes/groups.tsv'
    output:
        directory('data/A_representatives')
    container:
        'renv/renv.sif'
    script:
        '../scripts/A_extract-representatives.R'


rule A_checkm:
    input:
        "data/A_representatives/{genome}/genome.fna.gz"
    output:
        directory("data/A_checkm/{genome}")
    threads: 8
    container: 'checkm-genome:1.2.0--pyhdfd78af_0'
    shell:
        """
        #!/bin/bash
        # setup a tmp working dir
        tmp=$(mktemp -d)
        mkdir $tmp/ref
        cp {input} $tmp/ref/genome.fna.gz
        gunzip -c $tmp/ref/genome.fna.gz > $tmp/ref/genome.fna
        # run checkin
        checkm lineage_wf -t {threads} -x fna $tmp/ref $tmp/out > $tmp/stdout
        # prepare output folder
        mkdir -p {output}
        # copy results over
        cp -r $tmp/out/* {output}/
        cp $tmp/stdout {output}/checkm.txt
        # cleanup
        rm -rf $tmp
        """

def A_aggregate_genomes(wildcards,
                        x = 'data/A_representatives/{tax_bio}/genome.fna.gz):

    """
    List for all genomes the fasta sequence file (default, but adjustable by
    parameter x).
    """
    # make sure that the rule 'A_representatives' finished running
    chk = checkpoints.A_representatives.get().output[0]
    # load taxonomy table without dots in names
    df = pd.read_csv('data/A_representatives/taxonomy.tsv', sep = '\t')
    df.columns = df.columns.str.replace('.', '_')
    # list genome files
    return sample_wise(x, df)

rule A_chkm_summary:
    input:
        A_aggregate_genomes
    output:
        touch('data/A_checkm/checkm_summary.tsv')


