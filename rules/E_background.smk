
# Filter out gappy regions in alignemnt
# Then compute pairwise sequence ID and GC content for later stratification
# Problem: Output cannot contain a function -> checkpoint for E_clustal rule
checkpoint E_filter:
    input:
        #'data/D_search-seqs-aln/{region}.fna.gz'
        #lambda wild: D_aggregate(wild, 'data/D_search-seqs-aln/{region}.fna.gz')
        'data/D_search-seqs-aln/done.flag'
    output:
        #'data/E_search-filtered/{region}.fna.gz'
        directory('data/E_search-filtered/'),
        'data/E_search-gaps.png'
    container: 'renv/renv.sif'
    conda: 'renv'
    threads: 32
    script:
        '../scripts/E_filter.R'

def E_aggregate(wildcards, x):
    chk = checkpoints.E_filter.get().output
    # Check which regions where exported
    # (Could be few than in D_search-seqs due to filtering)
    xs = glob('data/E_search-filtered/*.fna.gz')
    xs = [ os.path.basename(i) for i in xs ]
    xs = [ i.split('.')[0] for i in xs ]
    # expand desired string
    df = pd.DataFrame({'region': xs})
    return sample_wise(x, df)

# convert to clustal format, surprisingly Biostrings seems not to do that
rule E_clustal:
    input:
        'data/E_search-filtered/{region}.fna.gz'
    output:
        'data/E_search-filtered/{region}.aln'
    conda: '../biopython.yml'
    shell:
        """
        # Needed in shell to work with conda
        python3 << EOF
from Bio import SeqIO
import gzip
with gzip.open('{input}', 'rt') as h:
    records = SeqIO.parse(h, 'fasta')
    count = SeqIO.write(records, '{output}', 'clustal')
    print('Converted %i records' % count)
EOF
        """

# Shuffle with di-nucleotide content
rule E_sissiz:
    input:
        'data/E_search-filtered/{region}.aln'
    output:
        'data/E_search-shuffled/{region}.txt'
    container:
        'sissiz/sissiz-3.0.sif'
    shell:
        """
        SISSIz --clustal -n 1        \
        --tstv --simulate            \
        --read_seeds=389650868,16063 \
        {input} > {output}     || true
        # the true and touch makes an empty file if SISSIz fails
        # (eg too large input)
        touch {output}
        """

rule E_collect:
    input:
        #lambda wild: E_aggregate(wild, 'data/E_search-shuffled/{region}.fna.gz')
        lambda wild: E_aggregate(wild, 'data/E_search-shuffled/{region}.txt')
    output:
        touch('data/E_search-shuffled/done.flag')

# Assess overall GC content, seq id, and dinucleotide frequencies
# before/after shuffeling
checkpoint E_stat:
    input:
        'data/E_search-shuffled/done.flag'
    output:
        'data/E_search-filtered-stat.tsv'
        #! side-effect: Convert SISSIz txt to
        #'data/E_search-shuffled/{region}.fna.gz'
    log: 'snakelogs/E_stat.txt'
    container: 'renv/renv.sif'
    conda: 'renv'
    threads: 32
    script:
        '../scripts/E_stat.R'

# custom aggregation of only those SISSIz outputs that were converted to fna.gz
def E_bg_models(wildcards, x):
    chk = checkpoints.E_stat.get().output
    # Check which regions where exported
    # (Could be few than in D_search-seqs due to filtering)
    xs = glob('data/E_search-shuffled/*.fna.gz')
    xs = [ os.path.basename(i) for i in xs ]
    xs = [ i.split('.')[0] for i in xs ]
    # expand desired string
    df = pd.DataFrame({'region': xs})
    return sample_wise(x, df)


