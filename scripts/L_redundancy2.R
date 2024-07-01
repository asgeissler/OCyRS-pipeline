# Part 2 (after RNAdistance run)
# Assess redundancy of the novel CRSs by
# 1. Phylogenetic similarity f contained species
# 2. Overlap in sequence overlaps (relative to genomic coordinates)
# 3. Alignment of the consensus structures

library(tidyverse)

# Colorblind-friendly palettes of the Color Universal Design
# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")

################################################################################
# continuously load intermediate results from part 1 of the script and make plots
################################################################################

dat.no.crs.regions <-
  'data/L1_redundancy/no.crs.regions.tsv' |>
  read_tsv()

p1.regionbar <-
  dat.no.crs.regions |>
  count(n) |>
  ggplot(aes(n, nn)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = nn),
            size = 5,
            nudge_y = 50) +
  scale_x_continuous(breaks = 0:max(dat.no.crs.regions$n)) +
  xlab('CRS in search region') +
  ylab('No. search regions') +
  theme_bw(18)
# p1.regionbar

################################################################################

potential.redundant <-
  'data/L1_redundancy/potential.redundant.tsv' |>
  read_tsv()
pr.jaccard <-
  'data/L1_redundancy/potential.redundant.jaccard.tsv' |>
  read_tsv()

p2.jaccard <-
  pr.jaccard |>
  ggplot(aes(jaccard)) +
  stat_ecdf() +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('Jaccard similarity of species shared\nbetween CRSs detected within\nthe same search region') +
  ylab('Empirical cumulative density') +
  theme_bw(18)
# p2.jaccard


################################################################################

rel.over <-
  'data/L1_redundancy/relative.overlaps.tsv' |>
  read_tsv()

p3.relover <-
  rel.over |>
  # ignore summetric A-B B-A comparisons for general overview of density
  filter(pos.x < pos.y) |>
  ggplot(aes(x)) +
  stat_ecdf() +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_vline(xintercept = .9,  color = cbbPalette[7]) +
  geom_hline(yintercept = .15, color = cbbPalette[7]) +
  xlab('Sequence overlap relative to\n shorter sequence length') +
  ylab('Empirical cumulative density') +
  theme_bw(18)

# p3.relover


################################################################################

aln.prop <-
  'data/L1_redundancy/alignment.proportion.tsv' |>
  read_tsv()

p4.alnprop <-
  aln.prop |>
  ggplot(aes(prop)) +
  stat_ecdf()  +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_vline(xintercept = .9, color = cbbPalette[6]) +
  geom_hline(yintercept = .6, color = cbbPalette[6]) +
  xlab('Proportion of alignment overlapped\n(with rel. overlap ≥ 0.9)') +
  ylab('Empirical cumulative density') +
  theme_bw(18)
# p4.alnprop

################################################################################

redundant.candidates <-
  'data/L1_redundancy/redundant.candidates.tsv' |>
  read_tsv()


p5.components <-
  redundant.candidates |>
  select(group, no.motifs) |>
  unique() |>
  count(no.motifs) |>
  ggplot(aes(no.motifs, n)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = n),
            size = 5,
            nudge_y = 2) +
  annotate(
    'label',
    3.5, 25,
    label = sprintf(
      'Total candidate groups: %g\nCandidate redundant motifs: %g',
      redundant.candidates |>
        select(group) |>
        unique() |>
        nrow(),
      redundant.candidates |>
        nrow()
    ),
    size = 5,
    hjust = 1
  ) +
  scale_x_continuous(breaks = 2:3) +
  xlab('CRS in candidate redundant group') +
  ylab('No. candidate redundant group') +
  theme_bw(18)
# p5.components

################################################################################

rc.seq <-
  'data/L1_redundancy/redundant.candidates.consensus.tsv' |>
  read_tsv()

test.pairs <-
  'data/L1_redundancy/RNAdistance.pairs.tsv'  |>
  read_tsv()

################################################################################
################################################################################
# Part 2: Processing RNAdistance resutls
################################################################################
################################################################################


rna.dist <-
  'data/L1_redundancy/RNAdistance.output.txt' |>
  read_lines() |>
  keep(str_detect, '^f:') |>
  str_remove('f: ') |>
  as.integer()

prelim <-
  test.pairs |>
  mutate(RNAdistance = rna.dist) |>
  arrange(RNAdistance) |>
  select(-c(consensus.x, consensus.y, keep.x, keep.y))

prelim |>
  write_tsv('data/L_rnadistance.tsv') 


################################################################################
# From manual inspection:
# - Cutoff 4, inclusive
# - Keep longer motif
# - Resolve ties by no.species

fdr <-
  'data/I_fdr.tsv' |>
  read_tsv()
  
motif.tax.pos <-
  'data/K_motif-tax-pos.tsv' |>
  read_tsv()
  
no.species <-
  motif.tax.pos |>
  select(motif, species) |>
  unique() |>
  count(motif, name = 'no.species')

dat <-
  prelim |>
  filter(RNAdistance <= 4) |>
  select(group, no.motifs, motif.x, motif.y) |>
  # mutate(row = 1:n()) |>
  # pivot_longer(c(motif.x, motif.y)) |>
  # count(value) |> count(n)
  # 3 motifs appear twice, fun with pairs ...
  left_join(fdr, c('motif.x' = 'motif')) |>
  left_join(fdr, c('motif.y' = 'motif')) |>
  left_join(no.species, c('motif.x' = 'motif')) |>
  left_join(no.species, c('motif.y' = 'motif'))

assertthat::are_equal(
  dat |>
    filter(alignment.len.x == alignment.len.y) |>
    filter(no.species.x == no.species.y) |>
    nrow(),
  0
)
# no remaining ties

redundant <-
  dat |>
  mutate(choose.x = (alignment.len.x > alignment.len.y) | (no.species.x > no.species.y)) |>
  transmute(
    superior          = ifelse(choose.x, motif.x, motif.y),
    superior.len      = ifelse(choose.x, alignment.len.x, alignment.len.y),
    superior.species  = ifelse(choose.x, no.species.x, no.species.y),
    # the to removed motif
    redundant         = ifelse(choose.x, motif.y, motif.x),
    redundant.len     = ifelse(choose.x, alignment.len.y, alignment.len.x),
    redundant.species = ifelse(choose.x, no.species.y, no.species.x)
  )

# redundant |>
#   select(superior, redundant) |>
#   tidygraph::tbl_graph(edges = _, directed = TRUE) |>
#   ggraph::ggraph() +
#   ggraph::geom_edge_link(arrow = arrow()) +
#   ggraph::geom_node_label(aes(label = name))
# 11 motifs to be removed

redundant |>
  write_tsv('data/L_redundant.tsv')

################################################################################

p6.rnadist <-
  prelim |>
  ggplot(aes(RNAdistance)) +
  stat_ecdf() +
  geom_vline(xintercept = 5, color = cbbPalette[4]) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('Comparison consensus structures\n(RNAdistance)') +
  ylab('Empirical cumulative density') +
  annotate(
    'label',
    120, .1,
    label = sprintf(
      'RNAdistance cutoff: %g\n(choose by manual data inspection)\nRedundant motifs: %g',
      4,
      redundant |>
        select(redundant) |>
        unique() |>
        nrow()
    ),
    size = 5,
    hjust = 1
  ) +
  theme_bw(18)
# p6.rnadist

################################################################################
# combine all plots

cowplot::plot_grid(
  p1.regionbar,
  p2.jaccard,
  p3.relover,
  p4.alnprop,
  p5.components,
  p6.rnadist,
  labels = 'AUTO',
  label_size = 20
)
ggsave('data/L_redundant.jpeg', width = 16, height = 11, dpi = 400)
