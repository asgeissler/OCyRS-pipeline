# Script purpose:
# Explore the data extracted from proGenomes and determine from shared
# eggNOG terms the orthology groups on which to anchor the structure search.

library(tidyverse)

# path.tax <- 'data/A_representatives/taxonomy.tsv'
# path.genes <- 'data/A_representatives/genes.tsv.gz'
# path.kegg <- 'data/A_representatives/kegg.tsv.gz'

path.tax <- unlist(snakemake@input[['tax']])
path.genes <- unlist(snakemake@input[['genes']])
path.kegg <- unlist(snakemake@input[['kegg']])


# path.out.fig <- 'data/B_OGs.jpeg'
# path.out.tbl <- 'data/B_OGs.tsv'

path.out.fig <- unlist(snakemake@output[['fig']])
path.out.tbl <- unlist(snakemake@output[['tbl']])

# make sure script output is placed in log file
log <- file(unlist(snakemake@log), open="wt")
sink(log)

################################################################################
# Loading
tax <- read_tsv(path.tax)
genes <- read_tsv(path.genes)
kegg <- read_tsv(path.kegg)

################################################################################
# Info about the representative genomes

print(paste(
 'Downloaded representative genomes:' ,
 nrow(tax)
))

print('Species with more than one genomes are:')

tax %>%
  count(species) %>%
  filter(n > 1) %>%
  arrange(desc(n)) %>%
  rename(multigenomes = n) %>%
  print

################################################################################
# overall distribution

tax %>%
  count(order) %>%
  arrange(desc(n)) %>%
  left_join(
    tax %>%
      select(order, species) %>%
      unique %>%
      count(order) %>%
      rename(species = n),
    'order'
  ) %>%
  mutate_at('order', str_remove, '^[0-9]+ ') %>%
  mutate_at('order', replace_na, 'unknown / unspecified') %>%
  rename(
    'Phylogenetic order' = order,
    'No. genomes' = n,
    'No. species' = species
  ) %>%
  print

################################################################################

################################################################################
# Overview stats on KEGG annotations

# total number of terms
print('From proGenomes available KEGG terms:')
kegg %>%
  select(db, term) %>%
  unique %>%
  count(db) %>%
  print()

# genes covered
print('The proportion of coding genes covered are')
kegg %>%
  select(db, tax.bio.gene) %>%
  unique %>%
  count(db) %>%
  mutate(
    covered.proporiton = 100 * n / (
      nrow(genes %>% filter(type == 'CDS'))
      )
  ) %>%
  print()


# average per gene
print('The average number of terms per gene and sd are:')
kegg %>%
  select(tax.bio.gene, db, term) %>%
  count(tax.bio.gene, db) -> foo
foo %>%
  group_by(db) %>%
  summarize(
    avg = mean(c(n, rep(0, nrow(genes) - length(n)))),
    sd = sd(c(n, rep(0, nrow(genes) - length(n))))
  ) %>%
  ungroup() %>%
  print()

print('The max number of terms per gene are')
foo %>%
  group_by(db) %>%
  summarize(max.n = max(n))


# singletons
print('The number singleton terms are')
kegg %>%
  count(db, term) -> foo
foo %>%
  filter(n == 1) %>%
  count(db) %>%
  rename(singletons = n) %>%
  left_join(
    foo %>%
      count(db) %>%
      rename(total = n),
    'db'
  ) %>%
  mutate(proporiton = singletons / total * 100)
  

# large terms
print('The largest terms aer')
foo %>%
  group_by(db) %>%
  slice_max(n)

################################################################################
################################################################################
# Inspect phylogeny distribution of order/family/genus/species
# Q: What is the largest sub-clade in which all genomes have at least one
#    copy of the orthologoue

# The genomes per clade
tax %>%
  select(tax.bio, order, family, genus, species) %>%
  gather('sub', 'txid', - tax.bio) %>%
  left_join(
    count(., sub, txid),
    c('sub', 'txid')
  ) %>%
  drop_na -> tax.sub


# check for overlaps
tax.sub %>%
  # filter(n > 3) %>%
  rename(genomes = n) %>%
  arrange(sub, txid) %>%
  left_join(
    kegg %>%
      select(- tax.bio.gene) %>%
      unique,
    'tax.bio'
  ) %>%
  drop_na() %>%
  count(db, term, title, sub, txid, genomes) %>%
  rename(term.genomes = n) -> tax.kegg.lap

    
# Only keep clades in which all genomes match
tax.kegg.lap %>%
  filter(genomes == term.genomes) %>%
  # largest clade size per term
  group_by(db, term, title) %>%
  summarize(max.clade = max(genomes)) %>%
  ungroup -> kegg.max.clade


################################################################################
# Determine unambigous terms form weighted relative gene count

kegg %>%
  filter(db == 'ko') -> ko

# weight genes by number of different terms they are part of
ko %>%
  count(tax.bio.gene) %>%
  rename(terms = n) %>%
  mutate(w = 1 / terms) %>%
  # sum up weights per term
  left_join(ko, 'tax.bio.gene') %>%
  group_by(term) %>%
  summarize(ws = sum(w)) %>%
  # compare to overall number of genes/geomes of term
  left_join(
    ko %>%
      group_by(term) %>%
      summarize(
        uniq.genes = length(tax.bio.gene),
        uniq.genomes = tax.bio %>%
          unique %>%
          length
      ),
    'term'
  ) %>%
  # add info on what the largest recalled clade was
  left_join(kegg.max.clade, 'term') %>%
  mutate_at('max.clade', replace_na, 0) -> dat.ko

dat.ko %>%
  ggplot(aes(uniq.genomes, uniq.genes)) +
  geom_point() +
  ylab('No. genes annotated in term') +
  xlab('No. of different genomes') +
  theme_bw(18) +
  ggtitle('Sizes of KEGG orthology terms') -> p1

dat.ko %>%
  ggplot(aes(max.clade)) +
  geom_histogram(bins = 60) +
  scale_x_continuous(breaks = seq(0, max(dat.ko$max.clade), 10)) +
  geom_vline(xintercept = c(5, 10, 30), color = 'red') +
  xlab('Largest sub-clade covered by term') +
  theme_bw(18) +
  ggtitle('Overlap of terms with taxonomic sub-clades',
          '(order, family, genus, species)') -> p2
  
dat.ko %>%
  mutate(
    max.c = cut(max.clade, c(-1, 0, 1, 4, 10, 30, Inf),
                include.lowest = TRUE) %>%
      fct_recode(
        'None' = "[-1,0]",
        '1' = "(0,1]",
        '2..4' = "(1,4]",
        '5..10' = "(4,10]",
        '11.30' = "(10,30]",
        '31+' = "(30,Inf]"
      )
  ) -> foo
foo %>%
  ggplot(aes(uniq.genomes, ws / uniq.genes,
             col = max.c)) +
  scale_color_manual(values = c(
      '#000000', #black
      "#D55E00", #organge
      "#56B4E9", # light blue
      "#0072B2", # dark blue
      "#009E73", # green
      "#CC79A7"  # pink
    ),
    name = 'Largest sub-clade covered by term'
  ) +
  geom_point(size = 2, alpha = 0.7) +
  xlab('No. of different genomes') +
  ylab('Weighted relative size') +
  theme_bw(18) +
  theme(legend.position = 'bottom') +
  ggtitle('KEGG orthology term uniqueness') -> p3

foo %>%
  ggplot(aes(max.c, ws / uniq.genes)) +
  geom_violin() +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_hline(yintercept = 0.85, color = 'blue') +
  xlab('Largest sub-clade covered by term') +
  ylab('Weighted relative size') +
  theme_bw(18) +
  ggtitle('Unambiguous orthologies') -> p4
  

cowplot::plot_grid(
  p1, p2, p3, p4,
  labels = 'AUTO',
  label_size = 18
)
ggsave(path.out.fig, width = 16, height = 10)

################################################################################
# Select as candidate ortholog groups (OGs) thar are 
# 1. Overlapping a clade of at least 5 genomes
# 2. Relative weight  size >.85

dat.ko %>%
  mutate(rel.size = ws / uniq.genes) %>%
  filter(
    rel.size >= .85,
    max.clade >= 5
  ) %>%
  select(term) %>%
  left_join(kegg, 'term') %>%
  select(term, title, tax.bio, tax.bio.gene) -> cand.og

################################################################################
# Sanity check uniqueness with pairwise overlaps

print('Overlaps between candidate OGs')
cand.og %>%
  count(term) %>%
  left_join(cand.og, 'term') %>%
  left_join(., ., c('tax.bio', 'tax.bio.gene')) %>%
  filter(term.x < term.y) %>%
  count(term.x, n.x, term.y, n.y) %>%
  rename(shared = n) %>%
  mutate(
    jac = shared / (n.x + n.y - shared)
  ) %>%
  pull(jac) %>%
  summary()


################################################################################
# repetition of overview stats for the candidate OGs

# overall no
cand.og %>%
  select(term) %>%
  unique %>%
  nrow -> i1


# genes covered
cand.og %>%
  select(tax.bio.gene) %>%
  unique %>%
  nrow %>%
  `*`(100) %>%
  `/`(genes %>% filter(type == 'CDS') %>% nrow) -> i2

# average per gene
cand.og %>%
  select(tax.bio.gene, term) %>%
  count(tax.bio.gene) %>%
  pull(n) -> foo

bar <- c(foo, rep(0, nrow(genes) - length(foo)))
paste(
  mean(bar),
  '±',
  sd(bar)
) -> i3 

# max terms per gene
i4 <- max(foo)

# affected genomes
cand.og %>%
  select(tax.bio) %>%
  unique() %>%
  nrow -> i5

# singletons, should be 0 but still make sure nothing went wrong
cand.og %>%
  count(term) %>%
  filter(n == 1) %>%
  nrow %>%
  assertthat::are_equal(0)
  
list(
  'No. candidate OGs' = i1,
  'Coding genes covered %' = i2,
  'per gene no. terms' = i3,
  'max terms per gene' = i4,
  '?uses all genomes' = i5 == nrow(tax)
) %>%
  paste(names(.), '=', .) %>%
  str_c(collapse = '\n') %>%
  cat()

# large terms
print('The sizes of the candidate OGs are:')
cand.og %>%
  count(term) %>%
  pull(n) %>%
  summary %>%
  print()

################################################################################
# Finally, save out

write_tsv(cand.og, path.out.tbl)

################################################################################
################################################################################
