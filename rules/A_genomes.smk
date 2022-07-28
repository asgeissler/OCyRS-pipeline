
rule A_download:
    output:
        'data/A_progenomes/genomes.fna.gz',
        'data/A_progenomes/genes.fna.gz',
        'data/A_progenomes/proteins.faa.gz',
        'data/A_progenomes/genes.tsv',
        'data/A_progenomes/groups.tsv',
        'data/A_progenomes/antibiotic.tsv'
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
        """

