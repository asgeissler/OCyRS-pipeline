
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
