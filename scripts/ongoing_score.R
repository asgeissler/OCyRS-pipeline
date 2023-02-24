# Load sequences of motif hits
# Check if these are part of the motif alignments
# Determine per motif score cutoff

library(tidyverse)
library(patchwork)

# library(corrplot)

# in.fdr <- 'data/I_fdr.tsv'
# in.scores <- 'data/H_scores.tsv'
# in.cmstat <- 'data/I_cmstat.tsv'
# in.rfam <- 'data/G_rfam-cmsearch.tsv.gz'

in.fdr <- unlist(snakemake@input[['fdr']])
in.scores <- unlist(snakemake@input[['scores']])
in.cmstat <- unlist(snakemake@input[['cmstat']])
in.rfam <- unlist(snakemake@input[['rfam']])


# out.cutoffs <- 'data/J_gathering-scores.tsv'

# Colorblind-friendly palettes of the Color Universal Design
# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2",
                "#FF0000",
                "#D55E00",
                "#CC79A7")

###############################################################################
# Load data

rfam <- read_tsv(in.rfam)
fdr <- read_tsv(in.fdr)
scores <- read_tsv(in.scores)
cmstat <- read_tsv(in.cmstat)

###############################################################################
# Build up motif categories identified by motif-dir pairs

bind_rows(
  scores %>%
    filter(str_starts(dir, 'E_search-shuffled_')) %>%
    select(motif, dir) %>%
    mutate(cat = 'Background motifs'),
  
  fdr %>%
    filter(RNAphylo.fdr > 10) %>%
    select(motif) %>%
    mutate(
      dir = 'D_search-seqs',
      cat = 'FDR > 10%'
    ),
  
  scores %>%
    filter(dir == 'D_search-seqs') %>%
    select(motif, dir) %>%
    mutate(cat = 'All biological motifs'),
  
  fdr %>%
    filter(RNAphylo.fdr <= 10) %>%
    select(motif) %>%
    mutate(
      dir = 'D_search-seqs',
      cat = 'FDR ≤ 10%'
    ),
  
  fdr %>%
    filter(RNAphylo.fdr <= 10) %>%
    select(motif) %>%
    semi_join(
      cmstat %>%
        filter(!str_starts(motif, 'RF')) %>%
        filter(bifurcations >= 1),
        c('motif')
    ) %>%
    mutate(
      dir = 'D_search-seqs',
      cat = 'FDR ≤ 10% & ≥ 1 bifurcation'
    ),
  
  fdr %>%
    filter(RNAphylo.fdr <= 10) %>%
    select(motif) %>%
    semi_join(
      scores %>%
        filter(dir == 'D_search-seqs') %>%
        filter(
          # Either covaryation or power
          ( expected / nbpairs * 100 > 10) |
            ( observed / nbpairs * 100 > 10)
        ),
        c('motif')
    ) %>%
    mutate(
      dir = 'D_search-seqs',
      cat = 'FDR ≤ 10% & either power or\ncovarying ≥ 10%'
    ),
  
  scores %>%
    filter(dir == 'G_rfam-bacteria-seeds') %>%
    select(motif, dir) %>%
    semi_join(
      rfam %>%
        select(motif = rf) %>%
        unique,
      'motif'
    ) %>%
    semi_join(
      scores %>%
        filter(dir == 'G_rfam-bacteria-seeds') %>%
        filter(
          # Either covaryation or power
          ( expected / nbpairs * 100 > 10) |
            ( observed / nbpairs * 100 > 10)
        ),
        c('motif')
    ) %>%
    mutate(cat = 'Cyanobacterial Rfam & either power or\ncovarying ≥ 10%'),
  
  scores %>%
    filter(dir == 'G_rfam-bacteria-seeds') %>%
    select(motif, dir) %>%
    semi_join(
      rfam %>%
        select(motif = rf) %>%
        unique,
      'motif'
    ) %>%
    mutate(cat = 'Cyanobacterial Rfam'),
  
  scores %>%
    filter(dir == 'G_rfam-bacteria-seeds') %>%
    select(motif, dir) %>%
    mutate(cat = 'Bacterial Rfam')
) -> cats


###############################################################################
# Build one table of scores with prettier names

scores %>%
  left_join(
    cmstat %>%
      select(motif, bifurcations, rel_entropy_cm, rel_entropy_hmm) %>%
      mutate(dir = ifelse(
        str_starts(motif, 'RF'),
        'G_rfam-bacteria-seeds',
        'D_search-seqs'
      )),
    c('motif', 'dir')
  ) %>%
  transmute(
    motif, dir,
    # CMfinder pipeline stats
    'Length, log10' = log10(alen + 1),
    'Sequences, log10' = log10(nseq + 1),
    'RNAphylo, log10' = log10(RNAphylo + 1),
    'hmmpair, log10' = log10(hmmpair + 1),
    # Ratio
    'log10(RNAphylo / Length)' = log10(RNAphylo / alen),
    # R-Scape
    'Avg. SI %' = avgid,
    'Paired positions %' = 2 * nbpairs / alen * 100,
    'Alignment power %' = expected / nbpairs * 100,
    'Covarying bps %' = observed / nbpairs * 100,
    # CMstat
    # consensus_residues_len, expected_max_hit_len,
    Bifurcations = bifurcations,
    'Rel. entropy CM' = rel_entropy_cm,
    'Rel. entropy hmm' =  rel_entropy_hmm
  ) -> score.dat


###############################################################################
# Bring togehter with categories for nices plots

cats %>%
  # preserve order for plotting
  mutate_at('cat', as_factor) %>%
  left_join(score.dat, c('dir', 'motif')) %>%
  gather('score', 'value', - c(dir, motif, cat)) %>%
  ggplot(aes(cat, value, fill = cat)) +
  geom_boxplot() +
  scale_fill_manual(
    values = cbbPalette,
    name = NULL
  ) +
  xlab(NULL) +
  ylab(NULL) +
  # scale_y_continuous(breaks = seq(0, 1, .1)) +
  facet_wrap(~ score, scales = 'free_y', ncol = 3) +
  theme_bw(16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'hide',
    axis.text.x = element_blank()
  ) -> p1

cats %>%
  mutate_at('cat', as_factor) %>%
  count(cat) %>%
  mutate(n2 = prettyNum(n, big.mark = ',')) %>%
  # preserve order for plotting
  ggplot(aes(cat, n, label = n2, fill = cat)) +
  geom_bar(stat = 'identity') +
  geom_text(nudge_y = 100) +
  scale_fill_manual(
    values = cbbPalette,
    name = NULL
  ) +
  xlab(NULL) +
  ylab('No motifs or Rfam families') +
  theme_bw(16) +
  guides(fill = guide_legend(nrow = 3, ncol = 3)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'bottom',
    axis.text.x = element_blank()
  ) -> p2

# p1 + p2
cowplot::plot_grid(p1, p2)
