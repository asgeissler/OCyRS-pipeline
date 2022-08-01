# This script extracts from a progenomes dump
# 1. the representative genomes
# 2. writes them out per assembly

library(tidyverse)
library(Biostrings)

library(conflicted)
conflict_prefer("filter", "dplyr")

# general prokaryotic info
path.rep <- unlist(snakemake@input[['rep']])
path.tax <- unlist(snakemake@input[['tax']])
# project specific data
path.genomes.seq <- unlist(snakemake@input[['genomes_seq']])
path.proteins.seq <- unlist(snakemake@input[['genes_seq']])
path.genes.seq <- unlist(snakemake@input[['genes_seq']])
path.genes <- unlist(snakemake@input[['genes']])
path.groups <- unlist(snakemake@input[['groups']])

# # general prokaryotic info
# path.rep <- 'data/A_progenomes/representatives.txt'
# path.tax <- 'data/A_progenomes/specI_lineageNCBI.tab'
# # project specific data
# path.genomes.seq <- 'data/A_progenomes/genomes.fna.gz'
# path.proteins.seq <- 'data/A_progenomes/proteins.faa.gz'
# path.genes.seq <- 'data/A_progenomes/genes.fna.gz'
# path.genes <- 'data/A_progenomes/genes.tsv'
# path.groups <- 'data/A_progenomes/groups.tsv'

# Output: Per representative genome 1 folder under
path.out <- unlist(snakemake@output)
# path.out <- 'data/A_representatives'
dir.create(path.out)

# Overview over URLs/APIs used in this script
url.eggnog <- 'http://eggnog5.embl.de/download/eggnog_5.0/e5.og_annotations.tsv'
url.kegg <- 'https://rest.kegg.jp/%s/%s'

################################################################################
# Load info on what are the (all) prokaryotic rep. genomes
xs <- read_lines(path.rep)

################################################################################
# load set of genomes this project is about

genomes <- readDNAStringSet(path.genomes.seq)
# memo: names are "tax.bioproject.chromosome"

tibble(
  tax.bio.chr = names(genomes),
  tax.bio = str_remove(tax.bio.chr, '\\.[^.]+$')
) %>%
  # only keep representative genomes
  filter(tax.bio %in% xs) -> genomes.dat

genomes.dat %>%
  mutate(path = sprintf('%s/%s/genome.fna.gz', path.out, tax.bio)) %>%
  group_by(path) %>%
  do(chrs = .$tax.bio.chr) -> foo

# Create subfolder per genome first
foo %>%
  pull(path) %>%
  dirname() %>%
  map(dir.create)

# save genome sequences
with(foo, map2(
  path, chrs,
  ~ writeXStringSet(genomes[.y], .x, compress = TRUE)
))

################################################################################
# repeat similarly for the gene and protein sequences

list(
  'genes' = readDNAStringSet(path.genes.seq),
  'proteins' = readAAStringSet(path.proteins.seq)
) -> dat

dat %>%
  map2(names(dat), function(x, y) {
    # x <- dat$genes ; y <- 'gene'
    # remove meta info in deflines
    names(x) <- str_remove(names(x), '\t.*$')
    # determine the corresponding genome
    tibble(
      x.tax.bio.gene = names(x),
      x.tax.bio = str_remove(x.tax.bio.gene, '\\.[^.]+$')
    ) %>%
      filter(x.tax.bio %in% genomes.dat$tax.bio) %>%
      mutate(path = sprintf('%s/%s/%s.fasta.gz', path.out, x.tax.bio, y)) %>%
      group_by(path) %>%
      do(towrite = .$x.tax.bio.gene) %>%
      with(map2(
        towrite, path,
        ~ writeXStringSet(x[.x], .y, compress = TRUE)
      ))
  })

################################################################################
# load meta information and tidy up a bit

genes <- read_tsv(path.genes) %>%
  mutate(tax.bio = str_remove(GENE_ID, '\\.[^.]+$')) %>%
  semi_join(genomes.dat, 'tax.bio')


# extract meta info
# scheme: key=list of words next_key=
# the str_extract below extracts the 'list of words' for the gene and product info
genes %>%
  select(tax.bio.gene = GENE_ID, ANNOTATION) %>%
  mutate(
    product = str_extract(
      ANNOTATION,
      '(?<=annotation product=)[^=]+(?= )'
    ),
    gene = str_extract(
      ANNOTATION,
      '(?<=gene=)[^=]+(?= )'
    )
  ) %>%
  select(- ANNOTATION) %>%
  # remove entries without either information
  filter(!is.na(product) | !is.na(gene)) -> meta

genes %>%
  select(
    tax.bio,
    tax.bio.gene = GENE_ID,
    tax.bio.chr = CONTIG_ID,
    type = TYPE,
    start = START,
    end = END,
    strand = STRAND,
    gc.content = `GC%`
  ) %>%
  left_join(meta, 'tax.bio.gene') -> genes2

# Write global coordinate list
write_tsv(genes2, sprintf('%s/genes.tsv.gz', path.out))

# And per genome
genes2 %>%
  mutate(path = sprintf('%s/%s/genes.tsv.gz', path.out, tax.bio)) %>%
  group_by(path) %>%
  do(i = list(select(., - tax.bio))) %>%
  ungroup() %>%
  with(map2(i, path, ~ write_tsv(.x[[1]], .y)))

################################################################################
# process orthology information

groups <- read_tsv(path.groups) %>%
  semi_join(genomes.dat, c('BIOSAMPLE_ID' = 'tax.bio'))

groups  %>%
  select(QUERY_NAME, EGGNOG_OGS) %>%
  # mutate(x = str_count(EGGNOG_OGS, ',')) %>%
  # pull(x) %>% summary()
  separate_rows(EGGNOG_OGS, sep = ',') %>%
  separate(EGGNOG_OGS, c('egg', 'taxid'), sep = '@') -> eggs

# Add human readable explanation texts
eggnog <- read_tsv(url.eggnog, col_names = c('taxid', 'egg', 'cog', 'text'))
# extracted from list on https://www.ncbi.nlm.nih.gov/research/cog/
tribble(
  ~ code, ~ category,
  'J', 'Translation, ribosomal structure and biogenesis',
  'A', 'RNA processing and modification',
  'K', 'Transcription',
  'L', 'Replication, recombination and repair',
  'B', 'Chromatin structure and dynamics',
  'D', 'Cell cycle control, cell division, chromosome partitioning',
  'Y', 'Nuclear structure',
  'V', 'Defense mechanisms',
  'T', 'Signal transduction mechanisms',
  'M', 'Cell wall/membrane/envelope biogenesis',
  'N', 'Cell motility',
  'Z', 'Cytoskeleton',
  'W', 'Extracellular structures',
  'U', 'Intracellular trafficking, secretion, and vesicular transport',
  'O', 'Posttranslational modification, protein turnover, chaperones',
  'X', 'Mobilome: prophages, transposons',
  'C', 'Energy production and conversion',
  'G', 'Carbohydrate transport and metabolism',
  'E', 'Amino acid transport and metabolism',
  'F', 'Nucleotide transport and metabolism',
  'H', 'Coenzyme transport and metabolism',
  'I', 'Lipid transport and metabolism',
  'P', 'Inorganic ion transport and metabolism',
  'Q', 'Secondary metabolites biosynthesis, transport and catabolism',
  'R', 'General function prediction only',
  'S', 'Function unknown'
) %>%
  arrange(code) -> cog.codes

eggs %>%
  mutate_at('taxid', as.numeric) %>%
  left_join(eggnog, c('egg', 'taxid')) %>%
  left_join(cog.codes, c('cog' = 'code')) %>%
  select(tax.bio.gene = QUERY_NAME,
         egg, taxlevel = taxid, egg.text = text,
         cog, cog.text = category) -> eggs2

# save as one large file for later gene matching
write_tsv(eggs2, sprintf('%s/eggnog.tsv.gz', path.out))

################################################################################
# Load KEGG pathway information in a readible form

# Collect terms
dat <- bind_rows(
  transmute(
    groups,
    tax.bio.gene = QUERY_NAME,
    db = 'pathway',
    term = KEGG_PATHWAY
  ) %>% drop_na,
  transmute(
    groups,
    tax.bio.gene = QUERY_NAME,
    db = 'module',
    term = KEGG_MODULE
  ) %>% drop_na,
  transmute(
    groups,
    tax.bio.gene = QUERY_NAME,
    db = 'ko',
    term = KEGG_KO
  ) %>% drop_na
) %>%
  separate_rows(term, sep = ',')
# look up descriptive text
dat %>%
  select(db) %>%
  unique() %>%
  mutate(q = sprintf(url.kegg, 'list', db)) %>%
  group_by(db) %>%
  do(i = read_tsv(.$q, col_names = c('term', 'title'))) %>%
  ungroup() %>%
  unnest(i) -> dat2

dat %>%
  mutate_at('term', str_remove, '^[a-z]+:') %>%
  inner_join(
    mutate_at(dat2, 'term', str_remove, '^[a-z]+:'),
    c('db', 'term')
  ) %>%
  mutate(tax.bio = str_remove(tax.bio.gene, '\\.[^.]+$')) -> dat3


# Save
write_tsv(dat3, sprintf('%s/kegg.tsv.gz', path.out))
dat3 %>%
  mutate(path = sprintf('%s/%s/kegg.tsv.gz', path.out, tax.bio)) %>%
  group_by(path) %>%
  do(i = list(select(., - tax.bio))) %>%
  ungroup() %>%
  with(map2(i, path, ~ write_tsv(.x[[1]], .y)))

################################################################################
# process taxonomy

read_tsv(
  path.tax,
  col_names = c('tax.bio', 'kingdom', 'phylum', 'class',
                'order', 'family', 'genus', 'species')
) %>%
  semi_join(genomes.dat, 'tax.bio') %>%
  write_tsv(sprintf('%s/taxonomy.tsv', path.out))


################################################################################
################################################################################