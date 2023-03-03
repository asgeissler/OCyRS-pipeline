# Purpose:
# Check for overlaps between motif homologs and cmsearch hits of
# (A) bacterial Rfam
# (B) RNIE terminator predictions
# Also, relate recall rates to motif score metriccs

library(tidyverse)
library(tidygraph)

in.rfam <- 'data/G_rfam-cmsearch.tsv.gz'
in.term <- 'data/G2_terminators.tsv.gz'
in.homologs <- 'data/J_motif-homologs.tsv'
in.cats = 'data/J_FDR-categories.tsv'

in.fdr <- 'data/I_fdr.tsv'
in.scores <- 'data/H_scores.tsv'

# Colorblind-friendly palettes of the Color Universal Design
# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")

################################################################################
# Load data

rfam <- read_tsv(in.rfam)
term <- read_tsv(in.term)
homologs <- read_tsv(in.homologs)
fdr <- read_tsv(in.fdr)
scores <- read_tsv(in.scores)


################################################################################
# as genomic range

list(
  'RNIE terminator' = term %>%
    select(- tax_bio, - gc),
  Rfam = rfam %>%
    select(- tax_bio, - gc) %>%
    rename(name = rf, rfam = name),
  'CMfinder motif' = homologs %>%
    rename(chr = tax.bio.chr, name = motif) %>%
    select(- c(tax_bio, gc, cms.row))
) %>%
  map2(names(.), ~ mutate(.x, type = .y)) %>%
  bind_rows() %>%
  mutate(row = 1:n()) %>%
  rename(seqnames = chr) %>%
  mutate(
    tmp = start,
    start = ifelse(start < end, start, end),
    end = ifelse(tmp < end, end, tmp)
  ) %>%
  select(- tmp) %>%
  plyranges::as_granges() %>%
  mutate(len = IRanges::width(.)) -> ranges

################################################################################
# check for overlaps all overlaps (also considering anti-sense)

ranges %>%
  mutate(strand2 = strand) %>%
  plyranges::join_overlap_intersect(., .) %>%
  # Exclude overlap with itself
  filter(row.x != row.y) %>%
  # Some helpful stats
  mutate(
    overlap = IRanges::width(.),
    jaccard = overlap / (len.x + len.y - overlap),
    x.rel = overlap / len.x,
    y.rel = overlap / len.y,
    orientation = ifelse(
      strand2.x == strand2.y,
      'sense',
      'anti-sense'
    )
  ) %>%
  as_tibble() -> overlaps

################################################################################
# First, filter out overlaps between rfam hits etc

overlaps %>%
  filter(row.x < row.y, type.x == type.y) -> over.type

cutoff.redundant <- 20

over.type %>%
  # ggplot(aes(jaccard, color = orientation)) +
  ggplot(aes(overlap, color = orientation)) +
  annotation_logticks(
    sides = 'b',
    short = unit(0.5, "cm"),
    mid = unit(0.8, "cm"),
    long = unit(1, "cm"),
  ) +
  scale_x_log10(minor_breaks = FALSE) +
  stat_ecdf(size = 1.2) +
  ggsci::scale_color_jama() +
  geom_vline(xintercept = cutoff.redundant, color = 'blue') +
  # geom_vline(xintercept = 40, color = 'red') +
  # scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('Overlap between cmsearch hits [bps]') +
  ylab('Emp. cum. density') +
  facet_wrap( ~ type.x) +
  theme_bw(18) +
  theme(legend.position = 'bottom')

ggsave('foo1.jpeg', width = 12, height = 6, dpi = 500)

# Determine overlaps above cutoff
over.type %>%
  filter(overlap >= cutoff.redundant) %>%
  select(from = row.x, to = row.y) %>%
  mutate_all(as.character) -> es
  
# Group of overlapping form graph
tbl_graph(edges = es, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(grp = group_components()) %>%
  as_tibble %>%
  mutate(row = as.integer(name)) %>%
  left_join(as_tibble(ranges), 'row') -> overlapping.hits

# Min e-value of overlapping group
overlapping.hits %>%
  group_by(grp) %>%
  summarize(min.e = min(evalue)) %>%
  ungroup %>%
  # In order to resolve ties, also choose max score
  left_join(overlapping.hits, 'grp') %>%
  filter(min.e == evalue) %>%
  group_by(grp, min.e) %>%
  summarize(max.score = max(score)) %>%
  ungroup %>%
  # Those with larger evalue than min per group are to be excluded
  left_join(overlapping.hits, 'grp') %>%
  filter( (min.e < evalue) | (score < max.score) ) %>%
  select(row) %>%
  unique -> to.exclude

ranges %>%
  filter(! (row %in% to.exclude$row)) -> ranges.nonredundant

overlaps %>%
  anti_join(to.exclude, c('row.x' = 'row')) %>%
  anti_join(to.exclude, c('row.y' = 'row')) -> overlaps.nonredundant

# length(ranges.nonredundant) / length(ranges)
# nrow(overlaps.nonredundant) / nrow(overlaps)

# Export non-redundant set
ranges.nonredundant %>%
  as_tibble() %>%
  rename(tax.bio.chr = seqnames) %>%
  select(- row, - len) %>%
  write_tsv('foo.tsv.gz')
  
################################################################################
# Overlap in non-redundant set ahead of recall

overlaps.nonredundant %>%
  filter(
    type.x == 'CMfinder motif',
    type.x != type.y
  ) -> over.between

cutoff.recall.overlap <- 5

over.between %>%
  ggplot(aes(overlap, color = orientation)) +
  stat_ecdf(size = 1.2) +
  ggsci::scale_color_jama() +
  # scale_x_continuous(breaks = c(seq(0, 50, 10), 100, 150)) +
  scale_x_log10(breaks = c(1, 5, 10, 50, 100, 150)) +
  annotation_logticks(sides = 'b') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_vline(xintercept = cutoff.recall.overlap, color = 'blue') +
  xlab('Overlap with CMfinder motifs [bps]') +
  ylab('Emp. cum. density') +
  facet_wrap( ~ type.y) +
  theme_bw(18) +
  theme(legend.position = 'bottom')

ggsave('foo2.jpeg', width = 12, height = 6, dpi = 500)

################################################################################
# Ahead of recall, query Rfam family types for more informative plots

'https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/bacteria.csv' %>%
  read_csv() %>%
  select(Family, Name = `Rfam ID`, type = `RNA type`) -> rfam.types

keep.type <- c('Cis-reg', 'antisense', 'CRISPR', 'ribozyme', 'rRNA', 'tRNA', 'sRNA')

rfam.types %>%
  mutate(
    type2 = type %>%
      map(function(i) {
        keep.type %>%
          map(~ str_detect(i, .x)) %>%
          unlist %>%
          which %>%
          first -> i
        keep.type[i]
      }) %>%
      unlist
  ) %>%
  mutate_at('type2', replace_na, 'other') %>%
  rename(type.full = type, type = type2) -> rfam.types

################################################################################

ranges.nonredundant %>%
  as_tibble() %>%
  count(name, type, name = 'no.seqs') -> no.seq


over.between %>%
  filter(overlap >= cutoff.recall.overlap) %>%
  left_join(rfam.types, c('name.y' = 'Family')) %>%
  # count(rfam.y) %>%
  # arrange(desc(n))
  mutate(
    type2 = case_when(
      type.y == 'RNIE terminator' ~ 'RNIE terminator',
      TRUE ~ paste('Rfam', type)
    )
  ) %>%
  select(motif = name.x, type = type.y, type2, name = name.y) %>%
  count(motif, type, type2, name, name = 'no.recall') %>%
  left_join(
    select(no.seq, motif = name, motif.homologs = no.seqs),
    'motif'
  ) %>%
  left_join(no.seq, c('name', 'type')) %>%
  mutate(
    recall = no.recall / no.seqs,
    precision = no.recall / motif.homologs,
    F1 = 2 * no.recall / ( 2 * no.recall + (no.seqs - no.recall) + (motif.homologs - no.recall))
  ) -> over.stat
over.stat %>%
  ggplot(aes(recall, precision, color = type2)) +
  # ggsci::scale_color_jco(name = NULL) +
  scale_color_manual(values = cbbPalette, name = NULL) +
  geom_point(size = 3, alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  xlab('recall (overlapping motifs rel. to\nno. RNIE/Rfam model hits)') +
  ylab('precision (overlapping motifs\nrel. to no. motif homologs') +
  annotation_logticks() +
  theme_bw(18) +
  theme(legend.position = 'bottom') -> p

p <- ggExtra::ggMarginal(p, type = 'histogram')
ggsave('foo3.jpeg', plot = p, width = 10, height = 10, dpi = 500)


###

ns <- sprintf('Overlapping RNA structure motifs (≥ %g bp)', cutoff.recall.overlap)
over.stat %>%
  select(motif, type2) %>%
  unique %>%
  count(type2) %>%
  arrange(desc(n)) %>%
  rename(!! ns:= n) %>%
  knitr::kable()

################################################################################
# Select the potential novel homologs


no.seq %>%
  filter(type == 'CMfinder motif') %>%
  select(- type) %>%
  rename(motif = name) %>%
  anti_join(over.stat, 'motif') -> potential.novel


ranges.nonredundant %>%
  filter(type == 'CMfinder motif') %>%
  filter(name %in% potential.novel$motif) %>%
  select(- type, - rfam) -> potential.novel.ranges

################################################################################
# Compare KEGG terms, associated pathways when comparing only the
# alignment sequences or also additional homologs

# Questions of interest:
# 1. Pathways associated via KO number
# 2. Pathways associated for genes +/- 500 bp away to include sequence homologs
# 3. When encluding cmsearch predicted homologs

# Q1
'https://rest.kegg.jp/link/pathway/ko' %>%
  read_tsv(col_names = c('ko', 'path')) %>%
  filter(str_detect(path, 'map')) %>%
  mutate_at('ko', str_remove, '^ko:') %>%
  mutate_at('path', str_remove, '^path:') -> ko.path

'https://rest.kegg.jp/list/pathway' %>%
  read_tsv(col_names = c('path', 'path_name')) %>%
  filter(str_detect(path, 'map')) %>%
  mutate_at('path', str_remove, '^path:') -> path.names

potential.novel %>%
  mutate(ko = str_remove(motif, '_.*$')) %>%
  left_join(ko.path, 'ko') %>%
  left_join(path.names, 'path') %>%
  select(path, path_name, motif) %>%
  unique %>%
  drop_na %>%
  count(
    path, path_name,
    name = 'Motifs associated via CMfinder search anchor'
  ) -> q1


# Q2

'data/A_representatives/kegg.tsv.gz' %>%
  read_tsv() -> gene.kegg

'data/A_representatives/genes.tsv.gz' %>%
  read_tsv() -> genes

# search in +/- 500
genes %>%
  select(
    seqnames = tax.bio.chr,
    start, end, strand,
    tax.bio.gene
  ) %>%
  plyranges::as_granges() %>%
  mutate(
    start = start - 500,
    end = end + 500
  ) -> genes500
  
# check overlaps
genes500 %>%
  plyranges::join_overlap_intersect(potential.novel.ranges) %>%
  as_tibble() -> genes500.novel.over


genes500.novel.over %>%
  left_join(gene.kegg, 'tax.bio.gene') %>%
  filter(db == 'pathway') %>%
  select(
    motif = name, tax.bio.gene, alignment.seq,
    path = term, path_name = title
  ) %>%
  unique -> genes500.motif.path
genes500.motif.path %>%
  filter(alignment.seq) %>%
  drop_na %>%
  count(
    path, path_name,
    name = 'Motifs with align. seq. near genes'
  ) -> q2

# Q3, repeated but without the alignment.seq filter
genes500.motif.path %>%
  drop_na %>%
  count(
    path, path_name,
    name = 'Motifs when incl. cmsearch homologs'
  ) -> q3


list(q1, q2, q3) %>%
  reduce(full_join, by = c('path', 'path_name')) %>%
  View

### Why does the number decrease for map00543?
potential.novel %>%
  mutate(ko = str_remove(motif, '_.*$')) %>%
  left_join(ko.path, 'ko') %>%
  left_join(path.names, 'path') %>%
  filter(path == 'map00543')

gene.kegg %>%
  filter(term == 'K19294') %>%
  select(tax.bio.gene) %>%
  left_join(gene.kegg) %>%
  filter(term == 'map00543')
# Problem: KEGG rest api lists KO-path associations not in proGenomes :*(

################################################################################

potential.novel %>%
  left_join(homologs, 'motif') %>%
  filter(str_detect(tax.bio.chr, '^32049')) %>%
  select(motif, tax.bio.chr, start, end, strand, alignment.seq)

potential.novel %>%
  left_join(homologs, 'motif') %>%
  filter(str_detect(tax.bio.chr, '^32049')) %>%
  select(motif) %>%
  left_join(fdr, 'motif')


################################################################################

