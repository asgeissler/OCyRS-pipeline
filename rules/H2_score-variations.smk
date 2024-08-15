# Test the various additional RNAphylo scoring options
# - For each reference tree
# - RNA structure alignment implied score

# Score the motifs with
# data/F_cmfinder/D_search-seqs/{region}/{region}.fna.motif.h1_1
# data/F_cmfinder/E_search-shuffled/{region}/{region}.fna.motif.h1_1
# with F_aggregate(wildcards, x)


rule H2_pscore:
    input:
        # here: {region} is split into {ortho}_{up/downstream side} to get the tree
        'data/H_symlink_{dir}/{ortho}_{side}.fna.motif.{number}.sto'
    output:
        'data/H2_scores/{dir}/{ortho}_{side}.fna.motif.{number}.txt'
    container:
        'cmfinder/cmfinder-0.4.1.18.sif'
    shell:
        """
        tmp=$(mktemp -d)
        echo "STO file:"  > {output}
        echo "{input}"   >> {output}
        
        # run for each reference tree
        for ref in reference-trees/*tree ; do
            # tree labels are: "txid.bioproject.name"
            # for RNAphylo, the '.name' has to be removed
            i=$(basename $ref .tree)
            sed -e 's,\.[^.:(),]*:,:,g' $ref > $tmp/$i.txt
            
            # compute RNAphylo per tree
            RNAPhylo --fragmentary --partition     \
                --partition-close-to-fuzzy-limit 3 \
                --suspicious-degen-nucs 2          \
                --ignore-all-gap                   \
                -t $tmp/$i.txt                     \
                {input}                            > $tmp/$i.pscore.txt
                
            # Store result in output
            echo "$i"                                    >> {output}
            grep 'Total .* posterior' $tmp/$i.pscore.txt >> {output}
        done

        # run once with tree implied by alignment
        perl $CMfinder/bin/ScoreMotif.pl {input} > $tmp/implied.txt
        
        echo "Implied tree"                       >> {output}
        grep 'RNAphylo score' $tmp/implied.txt    >> {output}

        rm -rf $tmp
        """



rule H2_all:
    input:
        lambda wild: F_aggregate(wild, 'data/H2_scores/{dir}/{filename}.txt'),
    output:
        touch('data/H2_scores/done.flag')


# this rule has to be run manually, since individual phylo scores just fail
# if not a single species overlap the reference tree.
# A clean Snakemake implementaiton is a bit too coplicates at this point
# use:
# bash run_local.sh H2_combine --allowed-rules H2_combine
rule H2_combine:
    input:
        #'data/H2_scores/done.flag',
        'data/H_scores.tsv'
    output:
        'data/H2_extra-scores.tsv',
        'data/H2_scores-comparisons.jpeg',
        'data/H2_shared-species.tsv',
        'data/H2_scores-species-tree.jpeg'
    container: 'renv/renv.sif'
    conda: 'renv'
    threads: 8
    script:
        '../scripts/H2_collect.R'
        
