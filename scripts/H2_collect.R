library(tidyverse)

library(furrr)

plan(multisession)
################################################################################

xs <-
  'data/H2_scores/*/*.txt' |>
  Sys.glob() 

################################################################################

extra.scores <-
  xs |>
  future_map(function(i) {
    # i <- 'data/H2_scores/D_search-seqs/K21929_upstream.fna.motif.h1_4.h2_3.h2_4.txt'
    i |>
      read_lines() |>
      discard(str_detect, 'Total pair posterior') |>
      str_remove('^RNAphylo score: ') |>
      str_remove('^Total RNA posterior ') |>
      str_remove('^data/H_symlink_') |>
      matrix(ncol = 2, byrow = TRUE) -> res
    set_names(res[, 2], res[, 1])
  }) |>
  bind_rows() |>
  separate(`STO file:`, c('dir', 'motif'), sep = '/') |>
  mutate_at(
    c('Chen_genomic', 'Chen_multigene', 'Cornet', 'Moore', 'Implied tree'),
    as.numeric
  )

extra.scores |>
  write_tsv('data/H2_extra-scores.tsv')

################################################################################

local.phylo <-
  'data/H_scores.tsv' |>
  read_tsv()

################################################################################

extra.scores |>
  pivot_longer(- c(dir, motif))  |>
  select(dir, motif, metric = name, rnaphylo = value) |>
  bind_rows(
    local.phylo |>
      transmute(dir, motif, metric = 'Local phylogeny', rnaphylo = RNAphylo) |>
      filter(dir != 'G_rfam-bacteria-seeds')
  ) |>
  filter(rnaphylo > 0) |>
  mutate(
    x = ifelse(dir == 'D_search-seqs', 'Biological motif', 'Random background'),
    metric2 = paste(x, metric, sep = ';') |>
      fct_reorder(rnaphylo),
    metric = paste(metric, sep = ';') |>
      fct_reorder(rnaphylo)
  ) -> foo

p1 <-
  foo |>
  ggplot(aes(rnaphylo, color = metric)) +
  stat_ecdf(size = 1.2) +
  scale_x_log10() +
  ggsci::scale_color_jama(name = 'Scoring relative to') +
  xlab('RNAphylo score\n(non-positive values omitted)') +
  ylab('Empirical cumulative density') +
  facet_wrap(~ x) +
  theme_bw(18)

p2 <-
  foo |>
  ggplot(aes(metric2, rnaphylo, color = metric, fill = metric)) +
  geom_violin() +
  geom_boxplot(fill = 'white', color = 'black', width = .1) +
  scale_y_log10() +
  ggsci::scale_color_jama(name = NULL) +
  ggsci::scale_fill_jama(name = NULL) +
  xlab('Scoring relative to\n(arranged by median)') +
  ylab('RNAphylo score\n(non-positive values omitted)') +
  scale_x_discrete(labels = foo |>
                     select(metric, metric2) |>
                     mutate_all(as.character) |>
                     unique() |>
                     with(set_names(metric, metric2))) +
  facet_wrap(~ x, scales = 'free_x') +
  theme_bw(18) +
  theme(legend.position = 'hide')


cowplot::plot_grid(p1, p2, ncol = 1)
ggsave('data/H2_scores-comparisons.jpeg', width = 22, height = 10)

################################################################################
# Check scores relatvie to species overlap

ref.spec <-
  'reference-trees/*.tree' |>
  Sys.glob() %>%
  set_names(., basename(.)) |>
  map(ape::read.tree) |>
  map(~ .$tip.label) |>
  map(str_remove, '\\..*$') |>
  map(unique)


extra.shared <-
  extra.scores |>
  select(dir, motif) |>
  mutate(
    region = str_remove(motif, '.fna.motif.*$'),
    path = 'data/F_cmfinder/%s/%s/%s' |>
      sprintf(dir, region, str_remove(motif, '\\.sto$'))
  ) |>
  with(future_pmap(list(dir, motif, path), function(dir, motif, path) {
    # path <- "data/F_cmfinder/D_search-seqs/K00008_downstream/K00008_downstream.fna.motif.h1_3"
    path.spec <-
      path |>
      read_lines() |>
      keep(str_starts, '[0-9]+\\.') |>
      str_remove(' +.*$') |>
      str_remove('\\..*$') |>
      unique()
    shared <-
      ref.spec |>
      map(~ intersect(.x, path.spec)) |>
      map(length) |>
      unlist()
    
    tibble(
      tree = shared |>
        names() |>
        str_remove('.tree$'),
      shared = shared,
      no.motif = length(path.spec),
      dir = dir,
      motif = motif
    )
  })) |>
  bind_rows()

extra.shared |>
  write_tsv('data/H2_shared-species.tsv')

################################################################################

extra.scores |>
  pivot_longer(- c(dir, motif))  |>
  inner_join(extra.shared, c('dir', 'motif', 'name' = 'tree')) |>
  mutate(
    x = ifelse(dir == 'D_search-seqs', 'Biological motif', 'Random background'),
  ) |>
  filter(value > 0) |>
  filter(x == 'Biological motif') |>
  ggplot(aes(shared / no.motif, value, color = name)) +
  geom_point() +
  xlab('Proportion of species in motif shared with tree') +
  ylab('RNAphylo score\n(non-positive values omitted)') +
  theme_bw(18) -> p
  # facet_wrap(~ x) -> p
             
p <- ggExtra::ggMarginal(p, groupColour = TRUE, groupFill = TRUE)
ggsave('data/H2_scores-species-tree.jpeg', p, width = 12, height = 12)
