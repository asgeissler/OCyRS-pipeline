# check overlap between annotations of
# - ProGenomes
# - Rfam CMsearch
# - RefSeq
# for the genomes of
# 
# - Synechococcus elongatus PCC 7942
# - Synechococcus sp. WH8102
# - Prochlorococcus marinus SS120
# - Prochlorococcus marinus MED4

library(tidyverse)
library(rentrez)
library(xml2)
library(plyranges)

################################################################################
# Semi-automatic ID lookup
tax <- read_tsv('data/A_representatives/taxonomy.tsv')

# match to tax.bio
xs <- c(
  # resolved by NCBI lookup of synonyms! Eg
  # https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=84588
  '84588.SAMEA3138327' = 'Synechococcus sp. WH8102',
  '1140.SAMN02598254' = 'Synechococcus elongatus PCC 7942',
  '167539.SAMN02603142' =  'Prochlorococcus marinus SS120',
  '59919.SAMEA3138209' =  'Prochlorococcus marinus MED4'
)

################################################################################
# Match chromosome IDs

proGenes <- read_tsv('data/A_representatives/genes.tsv.gz') %>%
  filter(tax.bio %in% names(xs))

proGenes %>%
  select(tax.bio, tax.bio.chr) %>%
  unique %>%
  mutate(chr = tax.bio.chr %>%
           strsplit('\\.') %>%
           map(3) %>%
           unlist) -> xs.chr

################################################################################
# Download Genbank Gene annotations

get.ref <- function(i) {
  qs <- entrez_search('nuccore', sprintf('%s [ACCN]', i))
  res <- entrez_fetch('nuccore', id = qs$ids, rettype = 'ft')
  
  res %>%
    read_tsv(
      skip = 1,
      col_names = c('start', 'end', 'info'),
      col_types = 'iic'
    ) %>%
    separate(info, c('key', 'value'), sep = '\t') %>%
    filter(is.na(value)) %>%
    mutate(
      strand = ifelse(start < end, '+',  '-'),
      tmp = start,
      start = ifelse(strand == '+', start, end),
      end = ifelse(strand == '+', end, tmp)
    ) %>%
    select(start, end, strand, type = key) %>%
    filter(type != 'gene') %>%
    mutate(chr = i)
}

xs.chr$chr %>%
  map(get.ref) %>%
  bind_rows() -> ref.genes

################################################################################
# Load rfam annotation

'data/G_rfam-cmsearch.tsv.gz' %>%
  read_tsv() %>%
  filter(tax_bio %in% names(xs)) %>%
  mutate(
    tmp = start,
    start = ifelse(strand == '+', start, end),
    end = ifelse(strand == '+', end, tmp)
  ) %>%
  select(- tmp) -> rfam

################################################################################
# Preapare comparison set

dat <- bind_rows(
  rfam %>%
    select(seqnames = chr, start, end, strand) %>%
    mutate(x = 'Rfam CMsearch'),
  proGenes %>%
    select(seqnames = tax.bio.chr, start, end, strand) %>%
    mutate(x = 'proGenomes annotation'),
  ref.genes %>%
    left_join(xs.chr, 'chr') %>%
    select(seqnames = tax.bio.chr, start, end, strand) %>%
    mutate(x = 'Genbank annotation')
)

################################################################################

dat %>%
  drop_na %>%
  # filter(!complete.cases(.))
  as_granges() %>%
  mutate(len = width(.), x.strand = strand) -> dat.ranges

dat.ranges %>%
  mutate(row = 1:n()) %>%
  join_overlap_intersect(., .) %>%
  mutate(
    m = ifelse(x.strand.x. == x.strand.y, 'sense', 'anti-sense'),
    jac = width / (len.x + len.y - width)
  ) -> dat2


################################################################################
dat2 %>%
  filter
################################################################################
################################################################################



