# Sensitive CMsearch of motifs against all ~500 bp long search regions
# 1. Default setting won't caputre the alignment sequences
# 2. Check for overlap between hits to determine some form of clustering
# 3. Ideally, clusters are correlated with a pathway

rule L_cmsearch:
    input:
        region = 'data/D_search-seqs/{region}.fna.gz',
        novel = 'data/J_novel/potentially-novel-motifs.tsv'
    output:
        directory('data/L_cmsearch-regions/{region}/')
    log: 'snakelogs/L_cmsearch-regions/{region}.txt'
    container: 'infernal\:1.1.4--pl5321hec16e2b_1'
    threads: 1
    shell:
        """
        # Query for the ~400 novel motifs (don't run all)
        cut -f 1 {input.novel}    | \
        tail -n +2                | \
        while IFS= read -r motif; do 
            # Select region and calibrated CM model 
            region=${{motif%%.fna.motif.*}}
            model=I_candidate-models/$region/$motif.cm
            # CMsearch with max sensetivity
            cmsearch --max --cpu {threads}    \
                --tblout {output}/$motif.txt  \
                $model            >> {log} 2>&1
        done
        """


rule L_cmsearch_agg:
    input:
        lambda wild: D_aggregate(wild, 'data/L_cmsearch-regions/{region}')
    output:
        touch('data/L_cmsearch-regions/done.txt')
        
        