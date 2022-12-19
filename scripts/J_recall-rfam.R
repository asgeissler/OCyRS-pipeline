# Purpose:
# Check for overlaps between motif homologs and Rfam hits.
# Check for recall of Rfam families also relative to motif scores.
# Use difference in score densities of recalling motifs and random
# background to select candidate motifs

library(tidyverse)

in.rfam <- 'data/G_rfam-cmsearch.tsv.gz'
in.motifs <- 'data/J_motif-homologs.tsv'
in.scores <- 'data/I_fdr.tsv'
in.all.scores <- 'data/H_scores.tsv'

################################################################################
# Load data
rfam <- read_tsv(in.rfam)
motifs <- read_tsv(in.motifs)
scores <- read_tsv(in.scores)
all.scores <- read_tsv(in.all.scores)


################################################################################
# as genomic range
rfam %>%
  select(
    seqnames = chr,
    start, end, strand,
    rf, name
  ) %>%
  mutate(
    tmp = start,
    start = ifelse(start < end, start, end),
    end = ifelse(tmp < end, end, tmp)
  ) %>%
  select(- tmp) %>%
  plyranges::as_granges() %>%
  mutate(rfam.len = IRanges::width(.)) -> rfam.range

motifs %>%
  select(
    seqnames = tax.bio.chr,
    start, end, strand,
    motif = name, score, evalue,
    motif.score.cutoff = min.seq.score
  ) %>%
  mutate(
    tmp = start,
    start = ifelse(start < end, start, end),
    end = ifelse(tmp < end, end, tmp)
  ) %>%
  select(- tmp) %>%
  plyranges::as_granges() %>%
  mutate(motif.len = IRanges::width(.)) -> motif.range

################################################################################
# check for overlaps

plyranges::join_overlap_intersect_directed(
  motif.range,
  rfam.range
) %>%
  mutate(
    shared = IRanges::width(.),
    jaccard = shared / (motif.len + rfam.len - shared),
    shared.rel = shared / motif.len
  ) -> over

over %>%
  as_tibble() %>%
  ggplot(aes(jaccard)) +
  stat_ecdf()
  ggplot(aes(shared.rel)) +
  stat_ecdf() +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_vline(xintercept = .5, color = 'red') +
  xlab('Rfam hit - motif cmsearch overlap\nrel. to motif hit length') +
  ylab('Emp. cum. density') +
  theme_bw(18)

over %>%
  as_tibble() %>%
  View
  nrow
  ggplot(aes(shared.rel, shared)) +
  geom_point(alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  annotation_logticks(sides = 'l') +
  scale_y_log10() +
  # scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_vline(xintercept = .5, color = 'red') +
  geom_hline(yintercept = 20, color = 'blue') +
  xlab('Rfam hit - motif cmsearch overlap\nrel. to motif hit length') +
  # ylab('Emp. cum. density') +
  theme_bw(18)
  
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
