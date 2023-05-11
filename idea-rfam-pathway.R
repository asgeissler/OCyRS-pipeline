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
  gather('k', 'v', - c( pathway, type )) %>%
  mutate(pathway = fct_reorder(pathway, v)) %>%
  ggplot(aes(pathway, v, fill = k)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~ type) -> p


ggsave('baz.jpeg', p)

################################################################################


