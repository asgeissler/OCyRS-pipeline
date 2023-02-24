# Load sequences of motif hits
# Check if these are part of the motif alignments
# Determine per motif score cutoff

library(tidyverse)
library(furrr)

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
# out.fig.eval <- 'data/J_eval.jpeg'
# out.fig.homologs <- 'data/J_motif-homologs.jpeg'

out.prop <- unlist(snakemake@output[['props']])
out.cutoffs <- unlist(snakemake@output[['cutoffs']])
out.homologs <- unlist(snakemake@output[['homologs']])

out.fig.jaccard <- unlist(snakemake@output[['fig_jaccard']])
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
  ylab('Empirical cumulative density') +
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
# Save out the alignment sequence recall, but do not filter here

dat %>%
  select(
    motif,
    'no.alignment.seqs' = no.seq,
    'no.cmsearch.found' = seq.recalled,
    'recalled.proportion' = prop
  ) %>%
  write_tsv(out.prop)

###############################################################################
# Ongoing for stefan question

dat %>%
  left_join(filter(scores, dir == 'D_search-seqs'), 'motif') %>%
  transmute(
    motif,
    prop = cut(
      prop,
      c(0, 50, 90, 100),
      include.lowest = TRUE
    ),
    # 'log10(RNAphylo)' = log10(RNAphylo),
    # 'log10(hmmpair)' = log10(hmmpair),
    'no. sequences' = nseq,
    'alignment length' = alen,
    'no. bps' = nbpairs,
    # 'significant co-variation' = observed
  ) %>%
  gather('score', 'value', - c(prop, motif)) %>%
  ggplot(aes(prop, value)) +
  geom_boxplot() +
  xlab('CMsearch recall of alignment sequences %') +
  facet_wrap(~ score, scales = 'free_y')

###############################################################################
# Determine E-value cutoff before gathering cutoff (massivly inflated)
cms %>%
  left_join(
    cmsearch.motif.overlap %>%
      filter(jacc == 1) %>%
      select(motif = motif.x, cms.row) %>%
      unique() %>%
      mutate(alignment.seq = TRUE),
    c('name' = 'motif', 'cms.row')
  ) %>%
  mutate_at('alignment.seq', replace_na, FALSE) %>%
  mutate(grp = ifelse(alignment.seq, 'Alignemnt seq.', 'Potential additional homolog')) %>%
  ggplot(aes(evalue, color = grp)) +
  stat_ecdf(size = 1.2) +
  ggsci::scale_color_jama(name = NULL) +
  # scale_x_log10(breaks = c(1e-60, 1e-30, 1e-10, 1e-6, 1e-2, 10)) +
  scale_x_log10() +
  xlab('CMsearch E-value for hits before filtering\n(log scaled axis)') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  # scale_y_log10() +
  geom_vline(xintercept = .1, color = 'black', linetype = 'dashed') +
  annotate('label', x = .1, y = 1, label = 'E=0.1') +
  # annotation_logticks(side = 'b') +
  ylab('Empirical cumulative density') +
  theme_bw(18) +
  theme(legend.position = 'bottom')

###############################################################################
# List potential homologs for all FDR10% motifs
# Use cmsearch screen to determine from minimal score of alignment sequence
# a sort of gathering score


cmsearch.motif.overlap %>%
  filter(jacc == 1) %>%
  group_by(motif = motif.x) %>%
  summarize(min.seq.score = min(score)) %>%
  ungroup -> cutoff

write_tsv(cutoff, out.cutoffs)

###############################################################################


cms %>%
  rename(motif = name) %>%
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
  # scale_x_log10(breaks = c(1e-30, 1e-20, 1e-10, 1e-6, 1e-2, 1)) +
  scale_x_log10() +
  xlab('CMsearch E-value for hits above motif gathering score\n(log scaled axis)') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_vline(xintercept = 1, color = 'black') +
  geom_vline(xintercept = c(0.01),
             color = 'black', linetype = 'dashed') +
  # annotation_logticks(side = 'b') +
  ylab('Empirical cumulative density') +
  theme_bw(18) +
  theme(legend.position = 'bottom')

ggsave(out.fig.eval, width = 12, height = 8)

###############################################################################
# List homologs

cand %>%
  filter(score >= min.seq.score) %>%
  filter(evalue <= .01) -> homologs

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
