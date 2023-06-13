library(tidyverse)

################################################################################

'data/J_novel/references_inside-of_intergenic_regions.tsv.gz' %>%
  read_tsv() -> rf

################################################################################
# Load curated KEGG annotation

kegg.path.ko <- 'https://rest.kegg.jp/link/pathway/ko' %>%
  read_tsv(col_names = c('term', 'path')) %>%
  mutate_all(str_remove, '^.*:')

kegg.paths <- 'https://rest.kegg.jp/list/pathway' %>%
  read_tsv(col_names = c('path', 'pathway'))

kegg.ko <- 'https://rest.kegg.jp/list/ko' %>%
  read_tsv(col_names = c('term', 'ortholog'))

################################################################################

rf.ko <- 
  rf %>%
  filter(type != 'RNIE terminator') %>%
  drop_na(region) %>%
  select(type, name, rfam, region) %>%
  unique %>%
  separate(region, c('ortholog', 'side'), sep = '_') %>%
  select(- side) %>%
  unique

################################################################################

rf.ko %>%
  left_join(kegg.path.ko, c('ortholog' = 'term')) %>%
  left_join(kegg.paths, 'path') %>%
  drop_na(pathway) -> rf.path

################################################################################

list(
  kegg.paths %>%
    left_join(kegg.path.ko, 'path') %>%
    count(pathway, name = 'KEGG orthologs in pathway'),
  
  rf.path %>%
    select(pathway, type, ortholog) %>%
    unique %>%
    count(pathway, type, name = 'KEGG orthologs in pathway with structure'),
  
  rf.path %>%
    select(type, name, pathway) %>%
    unique %>%
    count(pathway, type, name = 'Structures in pathway')
) %>%
  reduce(left_join) %>%
  # ignore pahtways without structrues
  drop_na(type) %>%
  mutate_if(is.numeric, replace_na, 0) -> dat

################################################################################

dat %>%
  mutate(prop = `KEGG orthologs in pathway with structure` / `KEGG orthologs in pathway` * 100) %>%
  select(pathway, type, prop) %>%
  spread(type, prop) %>%
  mutate_if(is.numeric, replace_na, 0) %>%
  ggplot(aes(Rfam, `Candidate motif`)) +
  geom_point() +
  ggtitle(
    'Orthologs with adjacent strucutres',
    'Proportion in number of KEGG orthologs annotated in pathway'
  ) +
  theme_bw(18)
  mutate_
  head

################################################################################

dat %>%
  gather('k', 'v', - c( pathway, type )) %>%
  mutate(pathway = fct_reorder(pathway, v)) %>%
  ggplot(aes(pathway, v, fill = k)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~ type) -> p


ggsave('baz.jpeg')

################################################################################
################################################################################
# OK. In addition I would like to see the number of search regions with
# (a) at least one cmsearch hit of the family, and
# (b) at least one cmsearch hit of the family and the specific KEGG pathway (count one for each KO search region).
# The log-odds ratio would be: log( b/a ).

rf %>%
  filter(type == 'Rfam') %>%
  mutate(term = str_remove(region, '_.*$')) %>%
  select(name, rfam, region, term) %>%
  unique %>%
  left_join(kegg.path.ko, 'term') %>%
  left_join(kegg.paths, 'path') %>%
  filter(str_starts(path, 'map')) -> rf.extra

list(
  
  rf.extra %>%
    select(rfam, name, region) %>%
    unique %>%
    count(rfam, name, name = 'A. Regions (KEGG orthologs but up/downstream specific) with structure'),
  
  rf.extra %>%
    select(pathway, rfam, name, region) %>%
    unique %>%
    count(pathway, rfam, name, name = 'B. Regions (KEGG orthologs but up/downstream specific) in pathway with structure')
  
) %>%
  reduce(left_join) %>%
  filter(
    `A. Regions (KEGG orthologs but up/downstream specific) with structure` >= 1,
    `B. Regions (KEGG orthologs but up/downstream specific) in pathway with structure` >= 1
  ) %>%
  mutate(
    log2.ratio.BA = log2(
      `B. Regions (KEGG orthologs but up/downstream specific) in pathway with structure` /
      `A. Regions (KEGG orthologs but up/downstream specific) with structure` 
    )
  ) %>%
  arrange(desc(log2.ratio.BA)) -> foo
foo %>%
  write_tsv('stefan-BA-ratio.tsv')

foo %>%
  ggplot(aes(log2.ratio.BA)) +
  stat_ecdf() +
  scale_y_continuous(breaks = seq(0, 1, .1))
