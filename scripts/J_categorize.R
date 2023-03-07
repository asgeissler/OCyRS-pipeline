# Inspect distribution of scores
# Further classify by variation and power

library(tidyverse)
# library(patchwork)
library(corrplot)

# in.fdr <- 'data/I_fdr.tsv'
# in.scores <- 'data/H_scores.tsv'
# in.cmstat <- 'data/I_cmstat.tsv'
# in.rfam <- 'data/G_rfam-cmsearch.tsv.gz'

in.fdr <- unlist(snakemake@input[['fdr']])
in.scores <- unlist(snakemake@input[['scores']])
in.cmstat <- unlist(snakemake@input[['cmstat']])
in.rfam <- unlist(snakemake@input[['rfam']])

# out.cor <- 'test.pdf'
# out.cut <- 'test1.png'
# out.dist <- 'test2.png'
# out.cats <- 'test2.tsv'

out.cor <- unlist(snakemake@output[['cor']])
out.cut <- unlist(snakemake@output[['cut']])
out.dist <- unlist(snakemake@output[['dist']])
out.cats <- unlist(snakemake@output[['cats']])

my.colors <- c(
  "Background motifs" = "#000000",
  "FDR > 10%" = "#999999",
  "All biological motifs" = "#A6CEE3", 
  "FDR ≤ 10%" = "#1F78B4",
  'High power\n(low covaryation)' = "#B15928",
  "Conserved sequence\n(low covaryation and power)" = "#FDBF6F",
  'High covaryation' = "#FF7F00",
  "Cyanobacterial Rfam" = "#6A3D9A",
  "Bacterial Rfam" = "#CAB2D6"
)
# RColorBrewer::brewer.pal(15, 'Paired')
# RColorBrewer::display.brewer.all()

###############################################################################
# Load data

rfam <- read_tsv(in.rfam)
fdr <- read_tsv(in.fdr)
scores <- read_tsv(in.scores)
cmstat <- read_tsv(in.cmstat)

###############################################################################
# Build up motif categories identified by motif-dir pairs

# Start collecting categories, but in parts, such to interject
# final catergories in "right place" for viz
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
    )
) -> cats1

  
bind_rows(
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
) -> cats2


###############################################################################
# Build one table of scores with prettier names

score.dat <- scores %>%
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
    # 'log10(RNAphylo / Length)' = log10(RNAphylo / alen),
    # R-Scape
    'Avgerage SI %' = avgid,
    'Paired positions %' = 2 * nbpairs / alen * 100,
    'Alignment power %' = expected / nbpairs * 100,
    'Covarying bps %' = observed / nbpairs * 100,
    # CMstat
    # consensus_residues_len, expected_max_hit_len,
    # Bifurcations = bifurcations,
    # 'Rel. entropy CM' = rel_entropy_cm,
    # 'Rel. entropy hmm' =  rel_entropy_hmm
  )

###############################################################################
# Determine covariaton and power cutoffs


bind_rows(cats1, cats2) %>%
  # preserve order for plotting
  mutate_at('cat', as_factor) %>%
  left_join(score.dat, c('dir', 'motif')) %>%
  select(
    motif, dir, cat,
    `Alignment power %`, `Covarying bps %`
  ) %>%
  gather('score', 'value', - c(dir, motif, cat)) %>%
  ggplot(aes(value, color = cat)) +
  stat_ecdf(size = 1.3) +
  scale_color_manual(
    values = my.colors[c(1, 2, 3, 4, 8, 9)],
    name = NULL
  ) +
  xlab(NULL) +
  ylab('Empirical cumulative density') +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  # scale_y_continuous(breaks = seq(0, 1, .1)) +
  facet_wrap(~ score, scales = 'free', strip.position = 'bottom') +
  geom_vline(xintercept = 20, color = 'red') +
  theme_bw(16) +
  theme(
    strip.placement = "outside",   # format to look like title
    strip.background = element_blank()
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank(),
    # legend.position = 'hide',
    # axis.text.x = element_blank()
  ) -> p.criteria

ggsave(filename = out.cut, plot = p.criteria,
       width = 16, height = 8,
       scale = 0.9,
       dpi = 400)

###############################################################################
# Create the 3 categories

fdr %>%
  filter(RNAphylo.fdr <= 10) %>%
  select(motif) %>%
  mutate(dir = 'D_search-seqs') %>%
  left_join(
    score.dat,
    c('motif', 'dir')
  ) -> fdr10

# Build iterativly, but do not in
fdr10 %>%
  mutate(cls = case_when(
    `Covarying bps %` >= 20 ~ 'High covaryation',
    `Alignment power %` >= 20 ~ 'High power\n(low covaryation)',
    TRUE ~ "Conserved sequence\n(low covaryation and power)"
  )) %>%
  select(motif, dir, cat = cls) %>%
  mutate_at(
    'cat', fct_relevel,
    'High power\n(low covaryation)',
    "Conserved sequence\n(low covaryation and power)",
    'High covaryation'
  ) %>%
  arrange(cat) %>%
  mutate_at('cat', as.character) -> cats3

###############################################################################
# Justification of selection due to correlation

cats <- bind_rows(cats1, cats3, cats2)

pdf(out.cor)
cats %>%
  filter(cat == 'FDR ≤ 10%') %>%
  left_join(score.dat, c('dir', 'motif')) %>%
  select(- c(motif, dir, cat))  %>%
  # GGally::ggpairs()
  cor %>%
  corrplot(
    method = 'square',
    order = 'AOE',
    type = 'lower',
    diag = FALSE,
    insig='blank',
    addCoef.col ='black',
    number.cex = 0.8
  )
dev.off()


###############################################################################
# Bring together with categories for nicer plots

cats %>%
  # preserve order for plotting
  mutate_at('cat', as_factor) %>%
  left_join(score.dat, c('dir', 'motif')) %>%
  gather('score', 'value', - c(dir, motif, cat)) %>%
  # ggplot(aes(str_remove(cat, '\\(.*\\)'), value, fill = cat)) +
  ggplot(aes(cat, value, fill = cat)) +
  geom_boxplot() +
  scale_fill_manual(
    values = my.colors,
    name = NULL
  ) +
  xlab(NULL) +
  ylab(NULL) +
  # scale_y_continuous(breaks = seq(0, 1, .1)) +
  facet_wrap(~ score, scales = 'free_y', ncol = 4) +
  # facet_wrap(~ score, scales = 'free_y') +
  theme_bw(16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'hide',
    axis.text.x = element_blank()
    # axis.text.x = element_text(angle = 90, hjust = 1)
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
    values = my.colors,
    name = NULL
  ) +
  xlab(NULL) +
  ylab('No. motifs or\nRfam families') +
  theme_bw(16) +
  # guides(fill = guide_legend(nrow = 3, size = 15)) +
  guides(fill = guide_legend(
    nrow = 2,
    label.hjust = 0.5, label.vjust = 1,
    label.position = 'bottom',
    byrow = TRUE,
    direction = 'vertical'
  )) +
  theme(
    # legend.box.spacing = unit(5, 'cm'),
    legend.key.width = unit(7, 'cm'),
    # legend.margin = margin(0, 0, 0, 0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'bottom',
    # legend.justification = c(0, 1),
    # axis.text.x = element_blank()
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) -> p2


cowplot::plot_grid(
  p1,
  p2,
  labels = 'AUTO',
  label_size = 18,
  rel_heights = c(1, 1.3),
  cols = 1
) -> p3


# p1 / p2 +
#   plot_annotation(tag_levels = 'A')

ggsave(filename = out.dist, plot = p3,
       width = 14, height = 12,
       bg = 'white',
       scale = 1.1,
       dpi = 400)

###############################################################################

cats3 %>%
  mutate_at('cat', str_replace, '\n', ' ') %>%
  left_join(score.dat, c('dir', 'motif')) %>%
  rename(category = cat) %>%
  select(- dir) %>%
  write_tsv(out.cats)
