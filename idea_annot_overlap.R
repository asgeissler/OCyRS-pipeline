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
  # filter(!complete.cases(.))
  drop_na %>%
  as_granges() -> dat.ranges

# dat.ranges %>%
#   filter(x != 'Genbank annotation') %>%
#   reduce_ranges_directed() %>%
#   mutate(x = 'proGenomes+Rfam') %>%
#   c(dat.ranges) %>%
dat.ranges %>%
  mutate(len = width(.), x.strand = strand) -> dat.ranges

dat.ranges %>%
  mutate(row = 1:n()) %>%
  join_overlap_intersect_directed(., .) %>%
  mutate(
    m = ifelse(x.strand.x == x.strand.y, 'sense', 'anti-sense'),
    jac = width / (len.x + len.y - width)
  ) -> dat2


################################################################################
# Compile overlap stats

dat2 %>%
  as_tibble() %>%
  filter(m == 'sense') %>%
  # Lookup species name
  left_join(
    xs.chr,
    c('seqnames' = 'tax.bio.chr')
  ) %>%
  left_join(
    tibble(
      tax.bio = names(xs),
      tax = xs
    ),
    'tax.bio'
  ) %>%
  # Ignore overlaps between annotations from same origin
  filter(x.x != x.y) %>%
  # Largest similarity of annotations in x to y
  group_by(tax, row.x, x.x, x.y, m) %>%
  summarize(jac = max(jac)) %>%
  ungroup %>%
  # de-duplicate pairs of same jaccard
  count(tax, x.x, x.y, m, jac) %>%
  # Compare to overall number of annotations in x
  left_join(
    dat.ranges %>%
      as_tibble() %>%
      # Lookup species name
      left_join(
        xs.chr,
        c('seqnames' = 'tax.bio.chr')
      ) %>%
      left_join(
        tibble(
          tax.bio = names(xs),
          tax = xs
        ),
        'tax.bio'
      ) %>%
      count(tax, x, name = 'nannot'),
    c('x.x' = 'x', 'tax')
  ) %>%
  #
  filter(jac >= .9) %>%
  group_by(tax, x.x, x.y, nannot) %>%
  summarize(n2 = sum(n)) %>%
  mutate(prop = n2 / nannot * 100) %>%
  # average
  group_by(x.x, x.y) %>%
  summarize(
    m = mean(prop),
    s = sd(prop)
  ) %>%
  ungroup %>%
  transmute(
    '↓ compared to →' = x.x,
    y = x.y,
    lab = sprintf('%.1f ±%.1f', m, s)
  ) %>%
  spread(y, lab, fill = '-')

################################################################################

