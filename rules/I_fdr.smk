# Aggregate scores and estimate FDRs
checkpoint I_fdr:
    input:
        'data/H_scores/done.flag',
        'data/H_scores.tsv'
    output:
        overview = 'data/I_overview.png',
        est = 'data/I_fdr.png',
        fdr = 'data/I_fdr.tsv',
        stat = 'data/I_alignment-stats.png'
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/I_fdr.R'


# Filter motifs by FDR 
def I_candidates_aggregate(wildcards, x):
    # make sure that the checkpoint finished
    chk = checkpoints.I_fdr.get().output
    # Load estimated FDR values and filter
    df = pd.read_csv('data/I_fdr.tsv', sep = '\t')
    mask = df['RNAphylo.fdr'] <= config['fdr_candidates']
    df2 = df.loc[mask]
    # expand wildecards
    return sample_wise(x, df2)


# Build covariance model for candidates
rule I_build:
    input:
        'data/F_cmfinder/D_search-seqs/{region}/{motif}'
    output:
        'data/I_candidate-models/{region}/{motif}.cm'
    container: 'infernal\:1.1.4--pl5321hec16e2b_1'
    threads: 4
    log: 'snakelogs/I_build/{region}/{motif}.txt'
    shell:
        """
        tmp=$(mktemp)
        rm $tmp # use file name, but prevent overwrite error of cmbuild
        cmbuild -n {wildcards.motif} -o {log} $tmp {input}
        cmcalibrate --cpu {threads} $tmp
        mv $tmp {output}
        """

# Compute stats for CM model
rule I_cmstat:
    input:
        'data/I_candidate-models/{region}/{motif}.cm'
    output:
        'data/I_cmstat/{region}/{motif}.txt'
    container: 'infernal\:1.1.4--pl5321hec16e2b_1'
    log: 'snakelogs/I_cmstat/{region}/{motif}.txt'
    shell:
        """
        cmstat {input} > {output}
        """

rule I_cmstat_agg:
    input:
        lambda wild: I_candidates_aggregate(wild, 'data/I_cmstat/{region}/{motif}.txt')
    output:
        'data/I_cmstat.tsv'
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
        '../scripts/I_cmstat.R'


rule I_build_cands:
    input:
        lambda wild: I_candidates_aggregate(wild, 'data/I_candidate-models/{region}/{motif}.cm')
    output:
        touch('data/I_candidate-models/done.cm')


# similar to G_cmsearch, but this time search for candidates
# Also, there is no gathering score to use
rule I_cmsearch:
    input:
        flag = 'data/I_candidate-models/done.cm',
        genome = 'data/A_representatives/{tax_bio}/genome.fna.gz'
    output:
        directory('data/I_cmsearch/{tax_bio}'),
        'data/I_cmsearch/{tax_bio}/run.done'
    log: 'snakelogs/I_rfam/{tax_bio}.txt'
    container: 'infernal\:1.1.4--pl5321hec16e2b_1'
    threads: 8
    shell:
        """
        mkdir -p {output[0]}
        for cm in data/I_candidate-models/*/*.cm ; do
            echo $cm >> {log} 2>&1
            p="{output[0]}/$(basename $cm .cm).txt"
            cmsearch --nohmmonly                \
                --tblout $p --cpu {threads}     \
                $cm {input.genome} >> {log} 2>&1
        done
        touch {output[1]}
        """


rule I_pergenome:
    input:
        lambda wild: A_aggregate_genomes(wild, 'data/I_cmsearch/{tax_bio}/run.done')
    output:
        touch('data/I_cmsearch/runs.done')


