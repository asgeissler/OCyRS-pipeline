# Inspect potential Pathway associations of motifs

library(tidyverse)

in.novel <- 'data/J_novel/potentially-novel-motifs.tsv'
in.pos <- 'data/J_novel/references_inside-of_intergenic_regions.tsv.gz'
in.motifs <- 'data/J_motif-aln-seq-pos.tsv'
in.cats <- 'data/J_FDR-categories.tsv'
in.tax <- 'data/A_representatives/taxonomy.tsv'

################################################################################
# Load curated KEGG annotation

kegg.path.ko <- 'https://rest.kegg.jp/link/pathway/ko' %>%
  read_tsv(col_names = c('term', 'path')) %>%
  mutate_all(str_remove, '^.*:')


kegg.paths <- 'https://rest.kegg.jp/list/pathway' %>%
  read_tsv(col_names = c('path', 'pathway'))
  
kegg.ko <- 'https://rest.kegg.jp/list/ko' %>%
  read_tsv(col_names = c('term', 'ortholog'))


ko.path <- kegg.ko %>%
  inner_join(kegg.path.ko, 'term') %>%
  inner_join(kegg.paths, 'path')

write_tsv(ko.path, 'data/K_ko-path.tsv')

################################################################################
# Load pipeline data

potential.novel <- read_tsv(in.novel)
motifs <- read_tsv(in.motifs)
cats <- read_tsv(in.cats)
tax <- read_tsv(in.tax)
pos <- read_tsv(in.pos)

################################################################################
# Association to pathways

potential.novel %>%
  mutate( motif.ko = str_replace(
    motif,
    '^(K[0-9]{5}_[updown]{2,4}stream).fna.motif.*$',
    '\\1'
  )) %>%
  separate(motif.ko, c('term', 'side'), sep = '_') %>%
  select(motif, category, term, side) %>%
  inner_join(ko.path, 'term') -> motif.path

write_tsv(motif.path, 'data/K_motif-path.tsv')
################################################################################
# Associaiton to positions and species

potential.novel %>%
  select(motif, category) %>%
  left_join(pos, c('motif' = 'name')) %>%
  select(- c(score, evalue, type, rfam, len, region)) %>%
  mutate(tax.bio = str_remove(seqnames, '\\.[^.]*$')) %>%
  left_join(tax, 'tax.bio') -> motif.tax.pos

write_tsv(motif.tax.pos, 'data/K_motif-tax-pos.tsv')

motif.tax.pos %>%
  select(- c(seqnames, start, end, width, strand)) %>%
  unique -> motif.tax

write_tsv(motif.tax, 'data/K_motif-tax.tsv')

################################################################################
# A small summary table

left_join(
  potential.novel %>%
    count(category, name = 'potential.novel.motifs'),
  motif.path %>%
    select(motif, category) %>%
    unique %>%
    count(category, name = 'motifs.with.pathway'),
  'category'
) -> foo

foo %>%
  mutate(category = 'In total') %>%
  group_by(category) %>%
  summarize_if(is.numeric, sum) %>%
  bind_rows(foo, .) %>%
  write_tsv('data/K_overview.tsv')

################################################################################