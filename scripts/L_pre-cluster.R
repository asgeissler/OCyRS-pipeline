# Sort of benchmark per region CMsearch and prepare for potential clustering

library(tidyverse)
library(plyranges)

################################################################################
# Load data

rel.cmsearch <-
  'data/L_cmsearch-rel-pos.tsv.gz' |>
  read_tsv()

rel.aln <-
  'data/L_motif-aln-rel-pos.tsv.gz' |>
  read_tsv()

################################################################################

rel.aln.ranges <-
  rel.aln |>
  dplyr::rename(seqnames = chr) |>
  as_granges() %>%
  mutate(
    w = width(.),
    row = 1:n()
  )

rel.cmsearch.ranges <-
  rel.cmsearch |>
  dplyr::rename(seqnames = chr) |>
  mutate(
    tmp = start,
    start = ifelse(strand == '+', start, end),
    end = ifelse(strand == '+', end, tmp)
  ) |>
  select(- tmp) |>
  as_granges() %>%
  mutate(
    w = width(.),
    row = 1:n()
  )

################################################################################
# Check overlap with alignment sequences and proportion of recall

rel.aln.ranges |>
  join_overlap_intersect_directed(rel.cmsearch.ranges) |>
  as_tibble() |>
  mutate(jaccard = width / (w.x + w.y - width)) -> aln.cmsearch.cmp

################################################################################
# Recall per Jaccard cutoff
c(1, .99, .9, .8) |>
  map(function(i){
    aln.cmsearch.cmp |>
      filter(motif.x == motif.y, region.x == region.y) |>
      filter(jaccard >= i) |>
      select(motif.x, row.x) |>
      unique() |>
      count(motif = motif.x, name = 'recalled.aln') |>
      mutate(cut = as.character(i))
  }) |>
  bind_rows() |>
  full_join(
    rel.aln |>
      count(motif, name = 'no.aln'),
    'motif'
  ) |>
  # filter(!complete.cases(.)) # good, all motifs have ≥1 alignment with perfect
  mutate(prop = recalled.aln / no.aln * 100) %>%
  ggplot(aes(prop, group = cut, color = cut)) +
  stat_ecdf(size = 1.2) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  ggsci::scale_color_jama(name = 'Jaccard similarity') +
  xlab('% Proportion of alignment sequences recalled') +
  ylab('Empirical cumulative density') +
  theme_bw(18) +
  theme(legend.position = 'bottom') -> p1



################################################################################
# Recall per E-value

c(10, 1, 0.05, 0.001) |>
  map(function(i){
    aln.cmsearch.cmp |>
      filter(motif.x == motif.y, region.x == region.y) |>
      filter(jaccard >= .9) |>
      filter(evalue <= i) |>
      select(motif.x, row.x) |>
      unique() |>
      count(motif = motif.x, name = 'recalled.aln') |>
      mutate(cut = as.character(i))
  }) |>
  bind_rows() |>
  full_join(
    rel.aln |>
      count(motif, name = 'no.aln'),
    'motif'
  ) |>
  # filter(!complete.cases(.)) # good, all motifs have ≥1 alignment with perfect
  mutate(prop = recalled.aln / no.aln * 100) %>%
  ggplot(aes(prop, group = cut, color = cut)) +
  stat_ecdf(size = 1.2) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  ggsci::scale_color_jama(name = 'E-value cutoff') +
  xlab('% Proportion of alignment sequences recalled\n(Jaccard similarity fixed ≥ 0.9)') +
  ylab('Empirical cumulative density') +
  theme_bw(18) +
  theme(legend.position = 'bottom') -> p2



################################################################################
# General Score~Evalue inspection

aln.cmsearch.cmp |>
  ggplot(aes(score, evalue)) +
  geom_hex(bins = 200) +
  geom_hline(
    yintercept = c(10, 1, .05, .001),
    color = ggsci::pal_jama()(4)
  ) +
  ggrepel::geom_label_repel(
    aes(x, y, label = txt),
    color = ggsci::pal_jama()(4),
    # direction = 'y',
    data = tribble(
      ~ x, ~ y, ~ txt,
      200, 10,    'E: 10',
      200, 1,     'E: 1',
      200, 0.05,  'E: 0.05',
      200, 0.001, 'E: 0.001'
    )
  ) +
  xlab('CMsearch score') +
  ylab('CMsearch E-value') +
  scale_y_log10() +
  theme_bw(18) -> p3


################################################################################
# Determine a sorts of gathering score

aln.cmsearch.cmp |>
  filter(motif.x == motif.y, region.x == region.y) |>
  filter(jaccard >= .9, evalue <= 0.001) |>
  group_by(motif = motif.x) |>
  summarize(
    gather.score = min(score),
    gather.e = max(evalue)
  ) |>
  ungroup() -> motif.gather



rel.cmsearch |>
  left_join(motif.gather, 'motif') |>
  filter(score >= gather.score) -> screen

screen |>
  mutate(motif.region = str_remove(motif, '\\.fna\\.motif.*$')) |>
  mutate(in.search = ifelse(region == motif.region,
                            'CMsearch.in', 'CMsearch.other')) |>
  count(motif, in.search) |>
  spread(in.search, n, fill = 0) |>
  mutate(
    total = CMsearch.in + CMsearch.other,
    # prop = CMsearch.other / total * 100,
    motif = fct_reorder(motif, total) |> fct_rev()
  ) |>
  pivot_longer(- motif) |>
  mutate(name = case_when(
    name == 'CMsearch.in' ~ 'Search region',
    name == 'CMsearch.other' ~ 'Other region',
    name == 'total' ~ 'Total hits'
  )) |>
  ggplot(aes(motif, value, color = name, group = name)) +
  geom_line(size = 1.2, alpha = 0.7) +
  xlab('Motifs (sorted by total no. hits)') +
  ylab('No. CMsearch hits\nrel. to "gathering score"') +
  ggsci::scale_color_jco(name = NULL) +
  theme_bw(18) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = 'bottom'
  ) -> p4
  


################################################################################
# Check for overlap between motifs

screen |>
  select(
    seqnames = chr, start, end, strand,
    motif
  ) |>
  mutate(
    tmp = start,
    start = ifelse(strand == '+', start, end),
    end = ifelse(strand == '+', end, tmp)
  ) |>
  select(- tmp) |>
  as_granges() %>%
  mutate(w = width(.)) -> screen.ranges

screen.ranges %>%
  # Not directed overlap !!!!!!!!
  join_overlap_intersect(., .) |>
  filter(motif.x < motif.y) |>
  as_tibble() |>
  mutate(jaccard = width / (w.x + w.y - width)) -> m2.lap

m2.lap |>
  ggplot(aes(jaccard)) +
  stat_ecdf(size = 1.2) +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('Jaccard similarity of overlap between\nCMsearch hits of two different motifs') +
  ylab('Empirical cumulative density') +
  theme_bw(18)

ggsave('data/L_motif-similarity.jpg')


m2.lap %>%
  filter(jaccard >= .9) |>
  count(motif.x, motif.y, name = 'shared') |>
  left_join(
    screen |>
      count(motif.x = motif, name = 'no.x'),
    'motif.x'
  ) %>%
  left_join(
    screen |>
      count(motif.y = motif, name = 'no.y'),
    'motif.y'
  ) |>
  mutate(jaccard = shared / (no.x + no.y - shared)) |>
  arrange(desc(jaccard)) |>
  write_tsv('data/L_motif-similarity.tsv')

################################################################################
################################################################################
# Export overlap figures as grid

cowplot::plot_grid(
  p1, p2, p3, p4,
  labels = 'AUTO', label_size = 16
)
ggsave('data/L_cmsearch-stat.jpg', 
       width = 16, height = 12)
