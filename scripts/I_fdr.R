library(tidyverse)

# Colorblind-friendly palettes of the Color Universal Design
# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")


# Input impplied by flags
in.stats <- Sys.glob('data/E_search-shuffled-stat*tsv')
in.scores <- 'data/H_scores.tsv'

# out.overview <- 'foo.png'
# out.est <- 'foo.png'
# out.fdr <- 'foo.tsv'
# out.stat <- 'foo.png'

out.overview <- unlist(snakemake@output[['overview']])
out.est <- unlist(snakemake@output[['est']])
out.fdr <- unlist(snakemake@output[['fdr']])
out.stat <- unlist(snakemake@output[['stat']])

################################################################################
# Prior scripts already checked that the GC/seq. id/di-nucleotide content
# did not change due to filtering/shuffeling, see
#  ../../OCyRS/OCyRS-pipeline/snakelogs/E_stat/*.txt
# assert that shuffling did not change properties (too much)

# Focus on GC/seq.id stats of the region for the stratification
in.stats[[1]] %>%
  read_tsv() %>%
  select(region, sid = avg.seqid.filtered, gc = gc.filtered) -> region.stat


################################################################################
# Bin per input region

xs <- c(.02, .2, .8, .98)
xs.sid <- c(.02, .2, .5, .8, .98)
# xs <- seq(0, 1, .01)

region.stat %>%
  pull(sid) %>%
  quantile(xs.sid) %>%
  round() -> sid.cut

region.stat %>%
  pull(gc) %>%
  # quantile(seq(0, 1, .01))
  quantile(xs) %>%
  round(2) -> gc.cut

region.stat %>%
  ggplot(aes(gc, sid)) +
  geom_point(alpha = .7) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  scale_x_continuous(breaks = seq(0, 1, .05)) +
  xlab('GC content') +
  ylab('Avg. seq. identity %') +
  # show grid
  geom_hline(yintercept = sid.cut, color = cbbPalette[3]) +
  geom_vline(xintercept = gc.cut,  color = cbbPalette[3]) +
  # label percent content row
  annotate(
    'text',
    x = .28,
    y = {
      off <- c(0, sid.cut)
      shift <- (c(sid.cut, 100) - off ) / 2
      off + shift
    },
    size = 5,
    label = sprintf('%.0f%%', {
      diff <- c(xs.sid, 1) - c(0, xs.sid)
      diff * 100
    }),
    color = cbbPalette[7]
  ) +
  # label percent content col
  annotate(
    'text',
    y = 10,
    x = {
      off <- c(min(region.stat$gc), gc.cut)
      shift <- (c(gc.cut, max(region.stat$gc)) - off ) / 2
      off + shift
    },
    size = 5,
    label = sprintf('%.0f%%', {
      diff <- c(xs, 1) - c(0, xs)
      diff * 100
    }),
    color = cbbPalette[8]
  ) +
  theme_bw(16) +
  ggtitle('Properties of search regions', 'after "gappy" column filtering') -> p

ggExtra::ggMarginal(p) -> p.regions
# p.regions

################################################################################
region.stat %>%
  transmute(
    region,
    gc = cut(gc, c(0, gc.cut, 1)),
    asid = cut(sid, c(0, sid.cut, 100))
  ) -> bins

################################################################################
bins %>%
  count(gc, asid) %>%
  mutate(w = ifelse(n < 100, 'white', 'black')) %>%
  ggplot(aes(gc, asid, fill = log(n), label = n)) +
  # ggplot(aes(gc, asid, fill = n, label = n)) +
  geom_tile() +
  geom_text(aes(color = I(w)), size = 5) +
  scale_fill_viridis_c() +
  xlab('GC content') +
  ylab('Avg. seq. identity %') +
  theme_bw(16) +
  ggtitle('Regions per bin') -> p.bin.region
# p.bin.region

################################################################################

scores <- in.scores %>%
  read_tsv() %>%
  mutate(data = case_when(
    dir == "D_search-seqs" ~ 'Biological sequences',
    dir == 'G_rfam-bacteria-seeds' ~ 'Rfam',
    dir == "E_search-shuffled_seed_007_42" ~ 'Shuffled#1',
    dir == "E_search-shuffled_seed_123_456" ~ 'Shuffled#2', 
    dir == "E_search-shuffled_seed_389650868_16063" ~ 'Shuffled#3',
    dir == "E_search-shuffled_seed_654_321" ~ 'Shuffled#4', 
    dir == "E_search-shuffled_seed_789_987" ~ 'Shuffled#5'
  )) %>%
  mutate(region = str_remove(motif, '.fna.motif.*$'))

################################################################################

scores %>%
  filter(data != 'Rfam') %>%
  select(region, data) %>%
  left_join(bins, 'region') %>%
  count(data, gc, asid) %>%
  filter(data == 'Biological sequences') %>%
  mutate(w = ifelse(n < 100, 'white', 'black')) %>%
  ggplot(aes(gc, asid, fill = log(n), label = n)) +
  geom_tile() +
  geom_text(aes(color = I(w)), size = 5) +
  scale_fill_viridis_c() +
  xlab('GC content') +
  ylab('Avg. seq. identity %') +
  theme_bw(16) +
  ggtitle('Detected motifs per bin') -> p.bin.motif
# p.bin.motif

scores %>%
  filter(data != 'Rfam') %>%
  select(region, data) %>%
  left_join(bins, 'region') %>%
  filter(data != 'Biological sequences') %>%
  count(gc, asid) %>%
  mutate(w = ifelse(n < 100, 'white', 'black')) %>%
  ggplot(aes(gc, asid, fill = log(n), label = n)) +
  geom_tile() +
  geom_text(aes(color = I(w)), size = 5) +
  scale_fill_viridis_c() +
  xlab('GC content') +
  ylab('Avg. seq. identity %') +
  theme_bw(16) +
  ggtitle('Random motifs per bin') -> p.bin.rand
# p.bin.rand

  

################################################################################
cowplot::plot_grid(
  p.regions,
  p.bin.region,
  p.bin.motif,
  p.bin.rand,
  labels = 'AUTO',
  label_size = 18
) -> p

ggsave(out.overview, plot = p,
       width = 14, height = 12,
       bg = 'white')

################################################################################
################################################################################

scores %>%
  filter(data != 'Rfam') %>%
  left_join(bins, 'region') %>%
  filter(RNAphylo > 0) %>%
  ggplot(aes(data, log(RNAphylo), color = data)) +
  geom_boxplot() +
  # ggsci::scale_color_igv(name = NULL) +
  scale_color_manual(values = rev(cbbPalette)[-1], name = NULL) +
  facet_grid(fct_rev(asid) ~ gc) +
  xlab(NULL) +
  ylab('log(RNAphylo), score < 0 omitted') +
  theme_bw(16) +
  theme(
    axis.text.x = element_blank(),
    legend.position = 'bottom'
  ) +
  ggtitle('RNAphylo scores per bin') -> p.bin.phylo

################################################################################
# Custom FDR estimation

# `ecdf` won't work because it will be under-estimating
my.fdr <- function(xs) {
  # xs is the random background
  step <- 1 / length(xs) * 100
  function(x) {
    # For a foreground observation x count no. background motifs with larger value
    up <- sum(xs > x)
    # The step size givew lower bound for FDR, each additional step increases
    # the error
    pmax(
      step,
      up * step
    )
  }
}


scores %>%
  filter(data != 'Rfam') %>%
  filter(data != 'Biological sequences') %>%
  left_join(bins, 'region') %>%
  filter(RNAphylo > 0) %>%
  group_by(gc, asid) %>%
  do(bg.model = my.fdr(.$RNAphylo)) -> bg.models


scores %>%
  filter(data == 'Biological sequences') %>%
  filter(RNAphylo > 0) %>%
  left_join(bins, 'region') %>%
  left_join(bg.models, c("gc", "asid")) %>%
  drop_na -> dat

with(
  dat,
  map2(
    bg.model, RNAphylo,
    function(fn, value) {
      fn(value)
    }
  )
) %>%
  unlist -> dat$fdr.est

################################################################################

dat %>%
  ggplot(aes(log10(RNAphylo), fdr.est)) +
  geom_point() +
  facet_grid(fct_rev(asid) ~ gc) +
  scale_y_log10(labels = scales::comma) +
  ylab('FDR estimate %') +
  theme_bw(16) -> p.est

dat %>%
  ggplot(aes(log10(RNAphylo), fdr.est)) +
  geom_point() +
  scale_y_log10(labels = scales::comma) +
  ylab('FDR estimate %') +
  theme_bw(16) -> p.est2


dat %>%
  ggplot(aes(fdr.est)) +
  stat_ecdf() +
  scale_x_log10(labels = scales::comma) +
  annotation_logticks(sides = 'b') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('FDR estimate %') +
  ylab('Cum. density') +
  theme_bw(16) -> p.est3

################################################################################
  
cowplot::plot_grid(
  p.bin.phylo,
  p.est,
  p.est2,
  p.est3,
  labels = 'AUTO',
  rel_heights = c(1.5, 1),
  label_size = 18
) -> p

ggsave(out.est, plot = p,
       width = 14, height = 12,
       bg = 'white')

################################################################################

dat %>%
  select(- c(bg.model, data, dir)) %>%
  transmute(
    motif,
    region,
    region.gc = gc, region.avgid = asid,
    RNAphylo, RNAphylo.fdr = fdr.est,
    hmmpair,
    no.seq = nseq, alignment.len = alen,
    avgid,
    no.bps = nbpairs,
    fraction.paired = nbpairs / alen * 100,
    expected.covary = expected,
    observed.covary = observed,
    alignment.power = expected / nbpairs * 100,
    fraction.covarying = observed / nbpairs * 100
  ) %>%
  write_tsv(out.fdr)

################################################################################

scores %>%
  filter(RNAphylo > 0) %>%
  transmute(
    data, motif,
    'RNAphylo, log10' = log10(RNAphylo),
    'hmmpair, log10' = log10(hmmpair),
    'Alignment length' = alen,
    'No. sequences in alignment' = nseq,
    'No. base-pairs' = nbpairs,
    'Fraction paired positions %' = nbpairs / alen * 100,
    'Average seq. id %' = avgid,
    'Alignment power %' = expected / nbpairs * 100,
    'Fraction covarying bps %' = observed / nbpairs * 100
  ) %>%
  mutate_at('data', str_remove, '#[0-9]$') %>%
  gather('k', 'value', - c(data, motif)) -> scores.concise
scores.concise %>%
  ggplot(aes(value, color = data)) +
  stat_ecdf(size = 1.5) +
  scale_color_manual(values = cbbPalette[-1], name = NULL) +
  xlab('Parameter value') +
  ylab('Emp. cum. density') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  facet_wrap(~ k, scales = 'free') +
  theme_bw(16) +
  theme(legend.position = 'bottom')

ggsave(out.stat, width = 12, height = 9)
