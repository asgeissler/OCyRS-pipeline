# Relative to RNIE models call terminators


rule G2_convert:
    input:
        'rnie-models/gene.cm'
    output:
        'data/G2_terminators/gene.cm'
    container: 'infernal\:1.1.4--pl5321hec16e2b_1'
    shell:
        """
        cmconvert -o {output} {input}
        """


rule G2_tmp:
    input:
        'data/A_representatives/{tax_bio}/genome.fna.gz'
    output:
        temp('data/tmp/{tax_bio}.fna')
    shell:
        """
        gunzip -c {input} > {output}
        """

rule G2_cmsearch:
    input:
        model = 'data/G2_terminators/gene.cm',
        genome = 'data/tmp/{tax_bio}.fna'
    output:
        'data/G2_terminators/{tax_bio}.txt'
    log: 'snakelogs/G2_cmsearch/{tax_bio}.txt'
    container: 'infernal\:1.1.4--pl5321hec16e2b_1'
    threads: 8
    shell:
        """
        #Legacy RNIE command for internal 0.81
        #cmsearch -T 14 -g --fil-no-qdb --fil-no-hmm --no-qdb --inside
        # In version 1.1.4 the correspondent command should be
        cmsearch -T 14 -g --nohmm       \
            --cpu 8                     \
            --tblout {output}           \
            {input.model} {input.genome} >> {log} 2>&1
        """

rule G2_pergenome:
    input:
        lambda wild: A_aggregate_genomes(wild, 'data/G2_terminators/{tax_bio}.txt')
    output:
        touch('data/G2_terminators/done.flag')

rule G2_combine:
    input:
        'data/G2_terminators/done.flag',
        'data/A_checkm/checkm_summary.tsv'
    output:
        tsv = protected('data/G2_terminators.tsv.gz'),
        fig = protected('data/G2_terminators.png')
    container: 'renv/renv.sif'
    conda: 'renv'
    threads: 8
    script:
        '../scripts/G2_combine_search.R'

