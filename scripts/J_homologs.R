# Load sequences of motif hits
# Check if these are part of the motif alignments
# Determine per motif score cutoff

library(tidyverse)
library(furrr)

library(corrplot)

# in.cmsearch <- 'data/J_cmsearch-collected.tsv'
# in.pos <- 'data/J_motif-aln-seq-pos.tsv'
# in.fdr <- 'data/I_fdr.tsv'
# in.scores <- 'data/H_scores.tsv'
# in.cmstat <- 'data/I_cmstat.tsv'
in.cmsearch <- unlist(snakemake@input[['cmsearch']])
in.pos <- unlist(snakemake@input[['pos']])
in.fdr <- unlist(snakemake@input[['fdr']])
in.scores <- unlist(snakemake@input[['scores']])
in.cmstat <- unlist(snakemake@input[['cmstat']])


# config.cutoff <- 10
config.cutoff <- unlist(snakemake@config[['fdr_candidates']])

# out.cutoffs <- 'data/J_gathering-scores.tsv'
# out.homologs <- 'data/J_motif-homologs.tsv'
# out.fig.jaccard <- 'data/J_overlaps.jpeg'
# out.fig.boxplot <- 'data/J_score-boxplots.jpeg'
# out.fig.cor <- 'data/J_score-cor.jpeg'
# out.fig.powcov <- 'data/J_high-power-covariation.jpeg'
# out.fig.eval <- 'data/J_eval.jpeg'
# out.fig.homologs <- 'data/J_motif-homologs.jpeg'

out.cutoffs <- unlist(snakemake@output[['cutoffs']])
out.homologs <- unlist(snakemake@output[['homologs']])

out.fig.jaccard <- unlist(snakemake@output[['fig_jaccard']])
out.fig.boxplot <- unlist(snakemake@output[['fig_boxplot']])
out.fig.cor <-  unlist(snakemake@output[['fig_cor']])
out.fig.powcov <- unlist(snakemake@output[['fig_powcov']])
out.fig.eval <-  unlist(snakemake@output[['fig_eval']])
out.fig.homologs <-  unlist(snakemake@output[['fig_homologs']])


cpus <- as.integer(unlist(snakemake@threads))
# cpus <- 20
plan(multisession, workers = cpus)

# Colorblind-friendly palettes of the Color Universal Design
# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")

###############################################################################
# Load coordinates

motif.seqs.pos <- read_tsv(in.pos)

cms <- read_tsv(in.cmsearch) %>%
  mutate(cms.row = 1:n())

# Split per genome
cms %>%
  rename(motif = name) %>%
  mutate(
    tmp = start,
    start = ifelse(start < end, start, end),
    end = ifelse(tmp < end, end, tmp)
  ) %>%
  select(- tmp) %>%
  group_by(tax_bio) %>%
  do(i = {
    grp <- .
    grp %>%
      select(seqnames = tax.bio.chr, start, end, strand,
             cms.row,
             motif, score, evalue) %>%
      plyranges::as_granges()
  }) %>%
  ungroup %>%
  with(set_names(i, tax_bio)) -> cms.ranges

###############################################################################
# No. sequences in motifs

motif.seqs.pos %>%
  count(motif, name = 'no.seq') -> no.seq

###############################################################################
# Check overlap between hits and alignment sequences

motif.seqs.pos %>%
  mutate(ms.row = 1:n()) -> ms.pos

names(cms.ranges) %>%
  future_map(function(i) {
    ms.pos %>%
      filter(str_starts(chr, i)) %>%
      rename(seqnames = chr) %>%
      plyranges::as_granges() %>%
      mutate(
        aln.len = IRanges::width(.),
        aln.strand = strand
      ) -> motif.range
    
    cms.ranges[[i]] %>%
      mutate(
        hit.len = IRanges::width(.),
        hist.strand = strand
      ) %>%
      # plyranges::join_overlap_intersect_directed(motif.range) %>%
      plyranges::join_overlap_intersect(motif.range) %>%
      as_tibble() %>%
      filter(motif.x == motif.y)
  }) %>%
  bind_rows() %>%
  mutate(jacc = width / (hit.len + aln.len - width)) -> cmsearch.motif.overlap

###############################################################################
# Check proportion of alignments recovered by CMsearch

cmsearch.motif.overlap %>%
  ggplot(aes(jacc)) +
  stat_ecdf() +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  xlab('Jaccard similarity cmsearch hit\nand motif alignment seq. pos.') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  ylab('Cum. emp. density') +
  theme_bw(18) -> p1

###############################################################################

cmsearch.motif.overlap %>%
  filter(jacc == 1) %>%
  select(
    motif = motif.x,
    cms.row, ms.row,
    aln.strand, hist.strand, 
  ) %>%
  filter(aln.strand == hist.strand) %>%
  select(motif, ms.row) %>%
  unique %>%
  count(motif, name = 'seq.recalled') %>%
  left_join(no.seq, 'motif') %>%
  mutate(prop = seq.recalled / no.seq * 100) -> dat

# Add excplicit 0 for motifs without a single recalled!
fdr <- read_tsv(in.fdr)
fdr %>%
  filter(RNAphylo.fdr <= config.cutoff) %>%
  select(motif) %>%
  unique %>%
  left_join(dat, 'motif') %>%
  # filter(!complete.cases(.)) %>%
  mutate_at('prop', replace_na, 0) %>%
  mutate(
    mf = motif %>%
      fct_reorder(prop) %>%
      fct_rev()) -> dat

# dat.prop.cutoff <- median(dat$prop)
dat.prop.cutoff <- 90

# plot curve of recall per motifs
dat %>%
  ggplot(aes(mf, prop, group = 1)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = dat.prop.cutoff, color = 'red') +
  geom_hline(yintercept = 50, color = 'blue') +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  xlab(sprintf(
    'Candidate motif (FDR <= %g%%)',
    config.cutoff
  )) +
  ylab('Poportion alignment\nseqeucnes recalled') +
  theme_bw(18) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) -> p2
# p2


dat %>%
  filter(prop >= dat.prop.cutoff) %>%
  select(motif) -> dat.better

print(sprintf(
  paste(
    'Of the n=%g motifs with FDR <%g%%, m=%g (%.1f%%) have a recall with',
    '(exact matches, same strand) above %.2g%%'
  ),
  nrow(dat),
  config.cutoff,
  nrow(dat.better),
  nrow(dat.better) / nrow(dat) * 100,
  dat.prop.cutoff
))
# [1] "Of the n=822 motifs with FDR <10%, m=462 (56.2%) have a recall with
#   (exact matches, same strand) above 90%"

cowplot::plot_grid(
  p1, p2,
  labels = 'AUTO',
  label_size = 18
)
ggsave(out.fig.jaccard,
       width = 16, height = 9,
       scale = 0.8)

###############################################################################
# Compare distribution of scores to substantiate selection by
# alignment sequence recall

scores <- read_tsv(in.scores)
cmstat <- read_tsv(in.cmstat)

scores %>%
  mutate(
    mode = case_when(
      dir == 'D_search-seqs' ~ 'CMfinder',
      dir == 'G_rfam-bacteria-seeds' ~ 'Rfam',
      TRUE ~ 'Background'
    )
  ) -> qual

bind_rows(
  qual %>%
    filter(mode != 'CMfinder') %>%
    # prevent mixing background motifs with cmstat
    mutate(motif = ifelse(
      mode == 'Background',
      paste0('bg', motif),
      motif
    )),
  qual %>%
    filter(mode == 'CMfinder') %>%
    anti_join(dat, 'motif') %>%
    mutate(mode = sprintf('FDR >%g%%', config.cutoff)),
  # qual %>%
  #   semi_join(dat, 'motif') %>%
  #   mutate(mode = sprintf('FDR ≤%g%%', config.cutoff)),
  qual %>%
    semi_join(
      dat %>%
        filter(prop >= dat.prop.cutoff) %>%
        select(motif),
      'motif'
    ) %>%
    mutate(mode = 'High aln. seq. recalling'),
  qual %>%
    semi_join(
      dat %>%
        filter(prop < dat.prop.cutoff) %>%
        filter(prop >= 50) %>%
        select(motif),
    ) %>%
    mutate(mode = 'Med. aln. seq. recalling'),
  qual %>%
    semi_join(
      dat %>%
        filter(prop < 50) %>%
        select(motif),
    ) %>%
    mutate(mode = 'Low aln. seq. recalling')
) %>%
  # pull(mode) %>% unique %>% dput
  mutate_at(
    'mode', fct_relevel,
    "Background",
    'FDR >10%',
    # 'CMfinder',
    # "FDR ≤10%",
    "Low aln. seq. recalling",
    "Med. aln. seq. recalling",
    "High aln. seq. recalling",
    "Rfam"
  ) %>%
  left_join(cmstat, 'motif') %>%
  transmute(
    mode, motif,
    # CMfinder pipeline stats
    'Length, log' = log10(alen + 1),
    'Sequences, log' = log10(nseq.x + 1),
    'RNAphylo, log' = log10(RNAphylo + 1),
    'hmmpair, log' = log10(hmmpair + 1),
    # Ratios
    # 'RNAphylo-Length, log-ratio' = `RNAphylo, log` / `Length, log`,
    # 'RNAphylo-Sequences, log-ratio' = `RNAphylo, log` / `Sequences, log`,
    # R-Scape
    'Avg. seq. id %' = avgid,
    'Fraction paired positions %' = 2 * nbpairs / alen * 100,
    'Alignment power %' = expected / nbpairs * 100,
    'Fraction covarying bps %' = observed / nbpairs * 100,
    # CMstat
    # consensus_residues_len, expected_max_hit_len,
    Bifurcations = bifurcations,
    'Rel. entropy CM' = rel_entropy_cm,
    'Rel. entropy hmm' =  rel_entropy_hmm
  ) -> xs


xs %>%
  gather('k', 'value', - c(mode, motif)) %>%
  ggplot(aes(mode, value, fill = mode)) +
  geom_boxplot() +
  scale_fill_manual(
    values = cbbPalette[c(1, 4, 3, 8, 7, 2)],
    # values = cbbPalette,
    name = NULL
  ) +
  xlab(NULL) +
  ylab('Parameter value') +
  # scale_y_continuous(breaks = seq(0, 1, .1)) +
  facet_wrap(~ k, scales = 'free_y', ncol = 3) +
  theme_bw(16) +
  theme(
    legend.position = 'bottom',
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

ggsave(out.fig.boxplot, width = 12, height = 16)

###############################################################################
# Show correlations between score metrics

xs %>%
  select(- c(motif, mode)) %>%
  drop_na() %>%
  cor -> xs.cor

# Order of score variables
# Do ahead of plot, such that the number coloration can be chosen
xs.ord <- corrplot::corrMatOrder(xs.cor, 'hclust')
xs.cor <- xs.cor[xs.ord, xs.ord]

jpeg(out.fig.cor, width = 3000, height = 3000, res = 400)
corrplot(
  xs.cor,
  order = 'hclust',
  method = 'square',
  type = 'upper',
  col = COL2('RdYlBu'),
  # addCoef.col = 'white',
  addCoef.col = xs.cor[upper.tri(xs.cor, diag = TRUE)] %>%
    abs %>%
    map( ~ ifelse(.x > .5, 'white', 'black')) %>%
    unlist
)
dev.off()

###############################################################################
# Argue via power cutoff for covariation cutoff

xs %>%
  # filter(mode == 'High aln. seq. recalling') %>%
  filter(`Alignment power %` >= 10) %>%
  # filter(Bifurcations >= 1) %>%
  ggplot(aes(`Fraction covarying bps %`, color = mode)) +
  stat_ecdf(size = 1.2) +
  scale_x_continuous(breaks = seq(0, 100, 10))  +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  ylab('Emp. cum. density') +
  scale_color_manual(
    values = cbbPalette[c(1, 4, 3, 8, 7, 2)],
    name = NULL
  ) +
  geom_vline(xintercept = c(25), linetype = 'dashed') +
  ggtitle('Min. alignement power 10%') +
  theme_bw(16) +
  theme(legend.position = 'bottom')
  
ggsave(out.fig.powcov, width = 8, height = 9)

###############################################################################
# Select candidate motifs

xs %>%
  filter(mode == 'High aln. seq. recalling') %>%
  filter(`Alignment power %` >= 10) %>%
  filter(`Fraction covarying bps %` >= 25) -> xs.cand

###############################################################################
###############################################################################
###############################################################################
# List potential homologs for candidates
# Use cmsearch screen to determine from minimal score of alignment sequence
# a sort of gathering score


cmsearch.motif.overlap %>%
  semi_join(xs.cand, c('motif.x' = 'motif')) %>%
  filter(jacc == 1) %>%
  group_by(motif = motif.x) %>%
  summarize(min.seq.score = min(score)) %>%
  ungroup -> cutoff

write_tsv(cutoff, out.cutoffs)

###############################################################################
# Determine E-value cutoff

cms %>%
  rename(motif = name) %>%
  semi_join(xs.cand, 'motif') %>%
  left_join(
    cmsearch.motif.overlap %>%
      filter(jacc == 1) %>%
      select(motif = motif.x, cms.row) %>%
      unique() %>%
      mutate(alignment.seq = TRUE),
    c('motif', 'cms.row')
  ) %>%
  mutate_at('alignment.seq', replace_na, FALSE) %>%
  left_join(cutoff, 'motif') %>%
  drop_na -> cand


cand %>%
  filter(score >= min.seq.score) %>%
  mutate(grp = ifelse(alignment.seq, 'Alignemnt seq.', 'Potential additional hit')) %>%
  ggplot(aes(evalue, color = grp)) +
  stat_ecdf(size = 1.2) +
  ggsci::scale_color_jama(name = NULL) +
  scale_x_log10(breaks = c(1e-30, 1e-20, 1e-10, 1e-6, 1e-2, 1)) +
  xlab('max. E-value alignment score per motif') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_vline(xintercept = 1, color = 'black') +
  geom_vline(xintercept = .01, color = 'black', linetype = 'dashed') +
  annotation_logticks(side = 'b') +
  ylab('Emp. cum. density') +
  theme_bw(18) +
  theme(legend.position = 'bottom')

ggsave(out.fig.eval, width = 12, height = 8)

###############################################################################
# List homologs

cand %>%
  filter(score > min.seq.score) %>%
  filter(evalue < .01) -> homologs

write_tsv(homologs, out.homologs)

###############################################################################
# A short figure on motifs per genome

homologs %>%
  select(motif, tax_bio) %>%
  unique %>%
  count(tax_bio) %>%
  ggplot(aes(n)) +
  stat_ecdf() +
  xlab('Motifs with at least one predicted homolog per genome') +
  scale_x_continuous(breaks = seq(5, 70, 5)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  ylab('Cum. emp. density') +
  theme_bw(18)

ggsave(out.fig.homologs, width = 9, height = 8)

###############################################################################