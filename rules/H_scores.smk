# Score the motifs with
# data/F_cmfinder/D_search-seqs/{region}/{region}.fna.motif.h1_1
# data/F_cmfinder/E_search-shuffled/{region}/{region}.fna.motif.h1_1
# with F_aggregate(wildcards, x)
#
# and the Rfam seeds in
# data/G_rfam-bacteria-seeds/{family}.sto
# with G_aggregate(wildcards, x)
#
# Challange: CMfinder scores will be relative to local phylogeny tree
# while Rfam is relative to the alignment tree.


# symlink CMfinder output to to simplify all rules below
rule H_symlink:
    input:
        #"data/F_cmfinder/D_search-seqs/KO123_upstream/KO123_upstream.fna.motif.h1_1"
        'data/F_cmfinder/{dir}/{region}/{region}.fna.motif.{number}'
    output:
        temp('data/H_symlink_{dir}/{region}.fna.motif.{number}.sto')
    shell:
        "ln -s {config[project_root]}/{input} {output}"
        

rule H_rscape:
    input:
        'data/{dir}/{file}.sto'
    output:
        directory('data/H_scores/rscape/{dir}/{file}/')
    container:
        'rscape/rscape-2.0.0.j.sif'
    shell:
        """
        mkdir -p {output}
        ln -s {config[project_root]}/{input} {output}/data.sto
        cd {output}
        R-scape -s data.sto > stdout.txt
        """


rule H_pscore:
    input:
        # here: {region} is split into {ortho}_{up/downstream side} to get the tree
        sto = 'data/H_symlink_{dir}/{ortho}_{side}.fna.motif.{number}.sto',
        tree = 'data/C_shrink/{ortho}/output.tree'
    output:
        'data/H_scores/cmscore/{dir}/{ortho}_{side}.fna.motif.{number}.txt'
    container:
        'cmfinder/cmfinder-0.4.1.18.sif'
    shell:
        """
        tmp=$(mktemp -d)
        RNAPhylo --fragmentary --partition     \
            --partition-close-to-fuzzy-limit 3 \
            --suspicious-degen-nucs 2          \
            --ignore-all-gap                   \
            -t {input.tree}                    \
            {input.sto}        > $tmp/pscore.txt

        hmmpair {input.sto}    \
            0.05 f 200 0 NULL  > $tmp/hmm.txt
        
        echo 'Pscore:'                             > {output}
        grep 'Total .* posterior' $tmp/pscore.txt >> {output}
        echo 'HMM:'                               >> {output}
        grep 'SCORE pairSumPf' $tmp/hmm.txt       >> {output}

        rm -rf $tmp
        """

# similar but for Rfam with the tree computed by the helper script
rule H_pscore2:
    input:
        'data/G_rfam-bacteria-seeds/{family}.sto'
    output:
        'data/H_scores/cmscore-Rfam/{family}.txt'
    container:
        'cmfinder/cmfinder-0.4.1.18.sif'
    shell:
        """
        perl $CMfinder/bin/ScoreMotif.pl {input} > {output}
        """


rule H_all:
    input:
        lambda wild: F_aggregate(wild, 'data/H_scores/rscape/H_symlink_{dir}/{filename}/'),
        lambda wild: G_aggregate(wild, 'data/H_scores/rscape/G_rfam-bacteria-seeds/{family}/'),
        lambda wild: F_aggregate(wild, 'data/H_scores/cmscore/{dir}/{filename}.txt'),
        lambda wild: G_aggregate(wild, 'data/H_scores/cmscore-Rfam/{family}.txt')
    output:
        touch('data/H_scores/done.flag')



rule H_combine:
    input:
        'data/H_scores/done.flag'
    output:
        'data/H_scores.tsv'
    container: 'renv/renv.sif'
    conda: 'renv'
    threads: 8
    script:
        '../scripts/H_collect_scores.R'
