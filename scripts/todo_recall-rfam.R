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

in.fdr <- 'data/I_fdr.tsv'
in.scores <- 'data/H_scores.tsv'

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
  geom_vline(xintercept = 20, color = 'blue') +
  # geom_vline(xintercept = 40, color = 'red') +
  # scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('Overlap between cmsearch hits [bps]') +
  ylab('Emp. cum. density') +
  facet_wrap( ~ type.x) +
  theme_bw(18) +
  theme(legend.position = 'bottom')

# Determine overlaps above cutoff
over.type %>%
  filter(overlap >= 20) %>%
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
  
################################################################################
# Overlap in non-redundant set ahead of recall

overlaps.nonredundant %>%
  filter(
    # row.x < row.y,
    type.x == 'CMfinder motif',
    type.x != type.y
  ) -> over.between
  # filter(row.x < row.y, type.x != type.y) -> over.between
  # filter(row.x < row.y) -> over.between

over.between %>%
  ggplot(aes(overlap, color = orientation)) +
  stat_ecdf(size = 1.2) +
  ggsci::scale_color_jama() +
  # scale_x_continuous(breaks = c(seq(0, 50, 10), 100, 150)) +
  scale_x_log10(breaks = c(1, 5, 10, 50, 100, 150)) +
  annotation_logticks(sides = 'b') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_vline(xintercept = 20, color = 'blue') +
  xlab('Overlap between cmsearch hits [bps]') +
  ylab('Emp. cum. density') +
  facet_wrap( ~ type.y) +
  theme_bw(18) +
  theme(legend.position = 'bottom')

################################################################################

ranges.nonredundant %>%
  as_tibble() %>%
  count(name, type, name = 'no.seqs') -> no.seq


over.between %>%
  filter(overlap >= 20) %>%
  count(rfam.y)
  mutate(
    type = case_when(
      type.y == 'RNIE terminator' ~ 'RNIE terminator',
      rfam.y
    )
  )
  View
  head
  select(motif = name.x, type, name = name.y) %>%
  head
  count(motif, type, name, name = 'no.recall') %>%
  left_join(
    select(no.seq, motif = name, motif.homologs = no.seqs),
    'motif'
  ) %>%
  left_join(no.seq, c('name', 'type')) %>%
  mutate(
    recall = no.recall / no.seqs,
    precision = no.recall / motif.homologs,
    F1 = 2 * no.recall / ( 2 * no.recall + (no.seqs - no.recall) + (motif.homologs - no.recall))
  ) %>%
  ggplot(aes(recall, precision)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks() +
  facet_wrap(~ type)

################################################################################

over %>%
  filter(shared.rel > .5) %>%
  as_tibble() %>%
  count(motif, rf, name, name = 'recalled') %>%
  left_join(
    motifs %>%
      count(motif = name, name = 'motif.pos'),
    'motif'
  ) %>%
  left_join(
    rfam %>%
      count(rf, name = 'rfam.pos'),
    'rf'
  ) -> recalls

recalls %>%
  # ggplot(aes(recalled / rfam.pos, recalled / motif.pos)) +
  # geom_point() +
  ggplot(aes( recalled / rfam.pos * 100)) +
  stat_ecdf() +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  ylab('Cum. emp. density') +
  xlab('% recall rfam hits') +
  geom_vline(xintercept = 5, color = 'red') +
  theme_bw(18) +
  ggtitle('Rfam-Motif recall', 'Motif hit overlaps > 50% of its length')

recalls %>%
  filter(recalled / rfam.pos * 100 > 5) %>%
  count(rf, name, name = 'No. motifs recalling > 5% of Rfam hits') %>%
  arrange(desc(`No. motifs recalling > 5% of Rfam hits`))

recalls %>%
  filter(recalled / rfam.pos * 100 > 5) %>%
  select(motif) %>%
  unique %>%
  nrow

################################################################################
scores %>%
  filter(RNAphylo.fdr <= 10) %>%
  transmute(
    motif,
    'log10 RNAphylo' = log10(RNAphylo),
    hmmpair,
    #RNAphylo.fdr, hmmpair, fraction.paired,
    alignment.power, fraction.covarying
  ) %>%
  left_join(
    recalls %>%
      filter(recalled / rfam.pos * 100 > 5) %>%
      select(motif) %>%
      unique %>%
      mutate(with.recall = TRUE),
    'motif'
  ) %>%
  mutate_at('with.recall', replace_na, FALSE) %>%
  gather('key', 'value', - c(motif, with.recall)) %>%
  ggplot(aes(value, color = with.recall)) +
  stat_ecdf() +
  ggsci::scale_color_jama() +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  ylab('Cum. emp. density') +
  xlab('score value') +
  facet_wrap(~ key, scales = 'free', ncol = 1) +
  theme_bw(18) +
  theme(legend.position = 'bottom')
  
################################################################################

################################################################################

all.scores %>%
  transmute(
    motif,
    dir = case_when(
      dir == 'D_search-seqs' ~ 'CMfinder',
      str_starts(dir, 'E_search-shuffled_seed') ~ 'Random CMfinder',
      dir == 'G_rfam-bacteria-seeds' ~ 'Rfam'
    ),
    logRNAphylo = log10(RNAphylo),
    log.hmmpair = log10(hmmpair),
    log.no.seq = log10(nseq),
    log.alignment.len = log10(alen),
    avgid,
    log.nbpairs = log10(nbpairs),
    fraction.paired = nbpairs / alen * 100,
    alignment.power = expected / nbpairs * 100,
    fraction.covarying = observed / nbpairs * 100
  ) -> foo
  
scores %>% 
  # TODO 10 from config
  filter(RNAphylo.fdr <= 10) %>%
  select(motif) %>%
  left_join(foo, 'motif') %>%
  mutate(dir = 'CMFinder, FDR≤10%') %>%
  bind_rows(foo) -> bar
bar %>%
  drop_na %>%
  select(- motif) %>%
  GGally::ggpairs(aes(color = dir))
