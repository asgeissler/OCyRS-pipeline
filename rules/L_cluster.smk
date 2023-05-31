# Sensitive CMsearch of motifs against all ~500 bp long search regions
# 1. Default setting won't caputre the alignment sequences
# 2. Check for overlap between hits to determine some form of clustering
# 3. Ideally, clusters are correlated with a pathway

rule L_cmsearch:
    input:
        region = 'data/D_search-seqs/{region}.fna.gz',
        novel = 'data/J_novel/potentially-novel-motifs.tsv'
    output:
        dir = directory('data/L_cmsearch-regions/{region}/'),
        flag = 'data/L_cmsearch-regions/{region}/done.flag'
    log: 'snakelogs/L_cmsearch-regions/{region}.txt'
    container: 'infernal\:1.1.4--pl5321hec16e2b_1'
    threads: 1
    shell:
        """
        # Query for the ~400 novel motifs (don't run all)
        cut -f 1 {input.novel}    | \
        tail -n +2                | \
        while IFS= read -r motif  ; do 
            # Select region and calibrated CM model 
            region=${{motif%%.fna.motif.*}}
            model=data/I_candidate-models/$region/$motif.cm

            echo "Motif $motif" >> {log}
            echo "Region $region" >> {log}
            echo "Model $model" >> {log}

            # CMsearch with max sensetivity
            cmsearch --max --cpu {threads}    \
                --tblout {output.dir}/$motif.txt  \
                $model {input.region}    >> {log} 2>&1
        done
        touch {output.flag}
        """

# Start CMsearch for novel motifs against search regions, but limited
# to only those regions that had motifs (~200 instead of >1k)
#  -> Less computationally intensive as D_aggregate
checkpoint L_helper:
    input:
        'data/J_novel/potentially-novel-motifs.tsv'
    output:
        'data/L_regions-todo.txt'
    shell:
        """
        cut -f 1 {input}           |
            tail -n +2             |
            sed 's,.fna.motif.*,,' |
            sort                   |
            uniq  > {output}
        """

def L_agg(w):
    # make sure that the checkpoint finished
    chk = checkpoints.L_helper.get().output
    # check which regions are needed
    with open('data/L_regions-todo.txt', 'r') as h:
        for reg in h:
            reg = reg.strip()
            yield 'data/L_cmsearch-regions/' + reg + '/done.flag'



rule L_cmsearch_agg:
    input:
        L_agg
    output:
        touch('data/L_cmsearch-regions/done.txt')



rule L_collect:
    input:
        'data/L_cmsearch-regions/done.txt',
        'data/L_regions-todo.txt',
        'data/J_novel/potentially-novel-motifs.tsv',
    output:
        'data/L_motif-aln-rel-pos.tsv.gz',
        'data/L_cmsearch-rel-pos.tsv.gz'
    container: 'renv/renv.sif'
    conda: 'renv'
    threads: 32
    script:
      '../scripts/L_rel-motif-pos.R'



rule L_precluster:
    input:
        'data/L_cmsearch-rel-pos.tsv.gz',
        'data/L_motif-aln-rel-pos.tsv.gz'
    output: 
        'data/L_motif-similarity.jpg',
        'data/L_motif-similarity.tsv',
        'data/L_cmsearch-stat.jpg', 
    container: 'renv/renv.sif'
    conda: 'renv'
    script:
      '../scripts/L_pre-cluster.R'
