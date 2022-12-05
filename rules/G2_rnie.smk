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


rule G2_cmsearch:
    input:
        model = 'data/G2_terminators/gene.cm',
        genome = 'data/A_representatives/{tax_bio}/genome.fna.gz'
    output:
        'data/G2_terminators/{tax_bio}.txt'
    container: 'infernal\:1.1.4--pl5321hec16e2b_1'
    threads: 8
    shell:
        """
        #Legacy RNIE command for internal 0.81
        #cmsearch -T 14 -g --fil-no-qdb --fil-no-hmm --no-qdb --inside
        # In version 1.1.4 the correspondent command should be
        cmsearch -T 14 -g --nohmm  \
            --cpu 8                \
            --tblout {output}      \
            {input.model} {input.genome}
        """

rule G2_pergenome:
    input:
        lambda wild: A_aggregate_genomes(wild, 'data/G2_terminators/{tax_bio}.txt')
    output:
        touch('data/G2_terminators/done.flag')

rule G2_combine:
    input:
        'data/G2_terminators/done.flag'
        'data/G_rfam-cmsearch/runs.done'
    output:
        protected('data/G2_terminators.tsv.gz')
    container: 'renv/renv.sif'
    conda: 'renv'
    threads: 8
    script:
        '../scripts/G_combine_search.R'

