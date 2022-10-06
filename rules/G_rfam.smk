
# Download all bacterial Rfam families
checkpoint G_download:
    output:
        directory('data/G_rfam-bacteria')
    params:
        rfamv = config['rfamv']
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/G_rfam-download.R'

checkpoint G_download_seeds:
    output:
        directory('data/C_rfam-bacteria-seeds')
    params:
        # The Rfam version to be used
        rfamv = config['rfamv']
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/G_rfam-download-seeds.R'


rule G_cmsearch:
    input:
        downflag = 'data/G_rfam-bacteria/download.done',
        genome = 'data/A_representatives/{tax_bio}/genome.fna.gz'
    output:
        directory('data/G_rfam-cmsearch/{tax_bio}')
    log: 'snakelogs/G_cmsearch-rfam/{tax_bio}.txt'
    container: 'infernal\:1.1.4--pl5321hec16e2b_1'
    threads: 8
    shell:
        """
        mkdir -p {output}
        for cm in data/G_rfam-bacteria/RF*.cm ; do
            echo $cm >> {log} 2>&1
            p="{output}/$(basename $cm .cm).txt"
            cmsearch --nohmmonly --rfam --cut_ga \
                --tblout $p --cpu 8              \
                $cm {input.genome} >> {log} 2>&1
        done
        touch {output}/run.done
        """

rule G_pergenome:
    input:
        lambda wild: A_aggregate_genomes(wild, 'data/G_rfam-cmsearch/{tax_bio}/run.done')
    output:
        touch('data/G_rfam-cmsearch/runs.done')

rule G_combine:
    input:
        'data/G_rfam-cmsearch/runs.done'
    output:
        protected('data/G_rfam-cmsearch.tsv.gz')
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/G_combine_search.R'

