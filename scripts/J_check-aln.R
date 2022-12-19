# Load sequences of motif hits
# Check if these are part of the motif alignments
# Determine per motif score cutoff

library(tidyverse)
library(furrr)


in.cmsearch <- 'data/J_cmsearch-collected.tsv'
in.pos <- 'data/J_motif-aln-seq-pos.tsv'
in.fdr <- 'data/I_fdr.tsv'
in.scores <- 'data/H_scores.tsv'
# in.cmsearch <- unlist(snakemake@input[['collected']])
# in.pos <- unlist(snakemake@input[['pos']])
# in.scores <- unlist(snakemake@input[['scores']])


config.cutoff <- 10
# config.cutoff <- unlist(snakemake@config[['fdr_candidates']])

# out.cutoffs <- 'data/J_score-cutoffs.tsv'
# out.homologs <- 'data/J_motif-homologs.tsv'
# out.fig <- 'data/J_overlaps.png'

# out.cutoffs <- unlist(snakemake@output[['cutoffs']])
# out.homologs <- unlist(snakemake@output[['homologs']])
# out.fig <- unlist(snakemake@output[['fig']])


# cpus <- as.integer(unlist(snakemake@threads))
cpus <- 20
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

# motif.seqs <- read_tsv(in.seq)
  
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
  mutate(prop = seq.recalled / no.seq * 100) %>%
  mutate(
    mf = motif %>%
      fct_reorder(prop) %>%
      fct_rev()) -> dat

dat.prop.cutoff <- median(dat$prop)

# plot curve of recall per motifs
dat %>%
  ggplot(aes(mf, prop, group = 1)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = dat.prop.cutoff, color = 'red') +
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


###############################################################################
# Compare distribution of scores to substantiate selection by data.prop.cutoff

scores <- read_tsv(in.scores)

scores %>%
  mutate(
    mode = case_when(
      dir == 'D_search-seqs' ~ 'CMfinder',
      dir == 'G_rfam-bacteria-seeds' ~ 'Rfam',
      TRUE ~ 'Background'
    )
  ) -> qual

bind_rows(
  qual,
  qual %>%
    semi_join(dat, 'motif') %>%
    mutate(mode = sprintf('FDR ≤%g%%', config.cutoff)),
  qual %>%
    semi_join(dat.better, 'motif') %>%
    mutate(mode = 'Aln. seq. recalling')
) %>%
  mutate_at(
    'mode', fct_relevel,
    "Background",
    'CMfinder',
    "FDR ≤10%",
    "Aln. seq. recalling",
    "Rfam"
  ) %>%
  transmute(
    mode, motif,
    'Length' = alen,
    'No. sequences' = nseq,
    # 'RNAphylo, log10' = log10(RNAphylo + 1),
    # 'hmmpair, log10' = log10(hmmpair + 1),
    # 'Fraction paired positions %, log10' = log10(nbpairs / alen * 100 + 1),
    # 'Alignment power %, log10' = log10(expected / nbpairs * 100 + 1),
    # 'Fraction covarying bps %, log10' = log10(observed / nbpairs * 100 + 1)
    'RNAphylo' = RNAphylo,
    'hmmpair' = hmmpair,
    'Fraction paired positions %' = nbpairs / alen * 100,
    'Alignment power %' = expected / nbpairs * 100,
    'Fraction covarying bps %' = observed / nbpairs * 100
  ) -> xs


xs %>%
  gather('k', 'value', - c(mode, motif)) %>%
  mutate_at('value', ~ log10(.x + 1)) %>%
  ggplot(aes(mode, value, color = mode)) +
  geom_violin() +
  # scale_y_log10() +
  facet_wrap(~ k, scales = 'free')

###############################################################################
xs2 %>%
  select(-motif) %>%
  mutate_at(
    'mode', fct_relevel,
    "Background",
    'CMfinder',
    "FDR ≤10%",
    "Aln. seq. recalling",
    "Rfam"
  ) %>%
  # mutate_if(is.numeric, ~ log10(.x + 1)) %>%
  mutate_at(
    c('Length', 'No. sequences', 'RNAphylo'),
    ~ log10(.x + 1)
  ) %>%
  GGally::ggpairs(aes(color = mode)) +
  scale_color_manual(
    values = cbbPalette[c(1, 4, 6, 7, 2)],
    name = NULL
  ) +
  scale_fill_manual(
    values = cbbPalette[c(1, 4, 6, 7, 2)],
    name = NULL
  )

# xs2 %>% 
xs %>% 
  select(mode, i = `Fraction paired positions %`) %>%
  group_by(mode) %>%
  do(is = list(.$i)) %>%
  ungroup %>%
  with(set_names(is, mode)) %>%
  map(1) -> foo

crossing(
 x =  names(foo),
 y =  names(foo)
) %>%
  filter(x < y) %>%
  head
t.test(foo$CMfinder, foo$Background, alternative = 'greater')
t.test(foo$`FDR ≤10%`, foo$CMfinder, alternative = 'greater')
t.test(foo$`Aln. seq. recalling`, foo$`FDR ≤10%`, alternative = 'greater')
t.test(foo$Rfam, foo$`Aln. seq. recalling`, alternative = 'greater')

t.test(foo$Rfam, foo$`Aln. seq. recalling`, alternative = 'two.sided')
###############################################################################

xs %>%
  mutate_at(
    'mode', fct_relevel,
    "Background",
    'CMfinder',
    "FDR ≤10%",
    "Aln. seq. recalling",
    "Rfam"
  ) %>%
  group_by(motif) %>%
  slice_max(mode, n = 1) %>%
  ungroup -> xs2

xs %>%
  filter(mode %in% c(
    "FDR ≤10%",
    "Aln. seq. recalling"
  )) %>%
  # filter(mode != 'Background') %>%
  # filter(mode != 'Rfam') %>%
  # mutate_if(is.numeric, ~ log10(.x + 1)) %>%
  mutate_at(
    c('Length', 'No. sequences', 'RNAphylo'),
    ~ log10(.x + 1)
  ) %>%
  drop_na %>%
  unite(
    'row',
    mode, motif, sep = ';'
  ) -> foo

foo %>%
  select(-row) %>%
  as.matrix %>%
  magrittr::set_rownames(foo$row) -> score.mat

score.mat %>%
  apply(2, scale) %>%
  # magrittr::set_rownames(rownames(score.mat)) -> score.scaled
  magrittr::set_rownames(rownames(score.mat)) %>%
  apply(1, scale) %>%
  t %>%
  magrittr::set_colnames(colnames(score.mat)) -> score.scaled

# boxplot(score.scaled)

score.pca <- prcomp(score.scaled)

autoplot(
  score.pca,
  colour = 'mode',
  alpha = 0.6,
  data = tibble(i = rownames(score.scaled)) %>%
    separate(i, c('mode', 'motif'), sep = ';') %>%
    mutate_at(
      'mode', fct_relevel,
      "Background",
      'CMfinder',
      "FDR ≤10%",
      "Aln. seq. recalling",
      "Rfam"
    ),
  # frame = TRUE, frame.type = 'norm',
  frame = TRUE, frame.type = 't',
  loadings = TRUE,
  loadings.colour = 'blue',
  loadings.label = TRUE,
  loadings.label.repel = TRUE,
  loadings.label.colour = 'black',
  loadings.label.fontface = 'bold',
  loadings.label.size = 7
) +
  scale_color_manual(
    values = cbbPalette[c(1, 4, 6, 7, 2)],
    name = NULL
  ) +
  scale_fill_manual(
    values = cbbPalette[c(1, 4, 6, 7, 2)],
    name = NULL
  )

###############################################################################
###############################################################################
  

xs %>%
  mutate_at(
    c('Length', 'No. sequences', 'RNAphylo'),
    ~ log10(.x + 1)
  ) %>%
  gather('k', 'value', - c(mode, motif)) %>%
  ggplot(aes(value, color = mode)) +
  stat_ecdf(size = 1.1) +
  scale_color_manual(
    values = cbbPalette[c(1, 4, 6, 7, 2)],
    name = NULL
  ) +
  xlab('Parameter value') +
  ylab('Emp. cum. density') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  facet_wrap(~ k, scales = 'free', ncol = 2) +
  theme_bw(16) +
  theme(legend.position = 'bottom')

xs %>%
  filter(mode == 'Aln. seq. recalling') %>%
  left_join(dat, 'motif') %>%
  View

xs %>%
  mutate(
    'Length, log10' = log10(Length + 1),
    'No. sequences, log10' = log10(`No. sequences` + 1),
    'RNAphylo, log10' = log10(RNAphylo + 1)
  ) %>%
  select(- c('Length', 'No. sequences', 'RNAphylo')) %>%
  select(- motif) %>%
  GGally::ggpairs(aes(color = mode)) +
  scale_color_manual(values = cbbPalette[c(1, 4, 6, 7, 2)]) +
  scale_fill_manual(values = cbbPalette[c(1, 4, 6, 7, 2)])

###############################################################################
###############################################################################
###############################################################################
###############################################################################
cmsearch.motif.overlap %>%
  filter(jacc > .99) %>%
  group_by(motif = motif.x) %>%
  summarize(min.seq.score = min(score)) %>%
  ungroup -> cutoff

write_tsv(cutoff, out.cutoffs)

###############################################################################

cms %>%
  rename(motif = name) %>%
  left_join(
    cmsearch.motif.overlap %>%
      filter(jacc > .99) %>%
      select(motif = motif.x, cms.row) %>%
      unique() %>%
      mutate(alignment.seq = TRUE),
    c('motif', 'cms.row')
  ) %>%
  mutate_at('alignment.seq', replace_na, FALSE) %>%
  left_join(cutoff, 'motif') %>%
  drop_na -> cand

cand %>%
  ggplot(aes(score / min.seq.score, - log10(evalue))) +
  geom_hex() +
  scale_fill_viridis_c() +
  xlab('ratio cmsearch hit score / \n max. motif alignment score') +
  theme_bw(18) -> p2

cmsearch.motif.overlap %>%
  filter(jacc > .99) %>%
  # pull(evalue) %>% summary
  ggplot(aes(evalue)) +
  geom_histogram(bins = 50) +
  scale_x_log10(breaks = c(1e-40, 1e-20, 1e-10, .05)) +
  xlab('max. E-value alignment score per motif') +
  theme_bw(18) -> p3

###############################################################################

cand %>%
  filter(score > min.seq.score) %>%
  filter(evalue < .05) -> homologs

write_tsv(homologs, out.homologs)

###############################################################################

full_join(
  homologs %>%
    count(motif = name, tax_bio) %>%
    group_by(motif) %>%
    summarize(avg.homo = mean(n)),
  motif.seqs %>%
    count(motif) %>%
    rename(no.seq = n),
  'motif'
) %>%
  ggplot(aes(avg.homo)) +
  stat_ecdf() +
  scale_x_log10() +
  xlab('Average no. predicted homolog motifs per genome\nafterwith E-value 0.05 and motif score') +
  annotation_logticks(sides = 'b') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  ylab('Cum. emp. density') +
  theme_bw(18) -> p4

###############################################################################

cowplot::plot_grid(
  p1, p2, p3, p4,
  labels = 'AUTO',
  label_size = 18
)

ggsave(out.fig, width = 20, height = 12)

###############################################################################