# Assess redundancy of the novel CRSs by
# 1. Phylogenetic similarity f contained species
# 2. Overlap in sequence overlaps (relative to genomic coordinates)
# 3. Alignment of the consensus structures

library(tidyverse)

library(plyranges)

library(conflicted)

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(dplyr::lag)

conflicted::conflicts_prefer(base::union)
conflicted::conflicts_prefer(base::setdiff)

# Colorblind-friendly palettes of the Color Universal Design
# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")

################################################################################
# Straight-forward pipeline data loading

motifs <-
  'data/K2_motifs.tsv' |>
  read_tsv() |>
  mutate(region = str_remove(motif, '.fna.*'))

pos <-
  'data/K_motif-tax-pos.tsv' |>
  read_tsv()

seqs <-
  'data/J_motif-aln-seqs.tsv.gz' |>
  read_tsv()

################################################################################
# Load exact positions of intergenic regions used for the screen

used.regions <-
  'data/D_search-seqs/*.fna.gz' |>
  Sys.glob() |>
  basename() |>
  str_remove('.fna.gz')

regions <-
  'data/J_novel/all_intergenic_regions.tsv.gz' |>
  read_tsv() |>
  # only keep entries that passed the min seq number criterion from
  filter(region %in% used.regions)

################################################################################
# Bar plot no motifs per search region

dat.no.crs.regions <-
  motifs |>
  count(region) |>
  right_join(tibble(region = used.regions), 'region') |>
  mutate_at('n', replace_na, 0)
  
p1.regionbar <-
  dat.no.crs.regions |>
  count(n) |>
  ggplot(aes(n, nn)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = nn),
            size = 5,
            nudge_y = 50) +
  xlab('CRS in search region') +
  ylab('No. search regions') +
  theme_bw(18)

################################################################################
################################################################################
# Part 1:
# Species similarity of motifs in the same search region
################################################################################
################################################################################

# motifs in region with min 2 CRS
potential.redundant <-
  dat.no.crs.regions |>
  filter(n >= 2) |>
  left_join(motifs, 'region') |>
  select(region, motif)

# the corresponding species per potential redundant motif (pr)
pr.species <-
  potential.redundant |>
  left_join(seqs, 'motif') |>
  select(- stockholm.seq) |>
  unique()

# the no. species per pr
pr.no.species <-
  pr.seq |>
  count(motif)

# all combinations of pairs, including self and symmetric case
pr.pairs <-
  potential.redundant %>%
  left_join(., ., 'region', relationship = 'many-to-many') |>
  left_join(pr.no.species, c('motif.x' = 'motif')) |>
  left_join(pr.no.species, c('motif.y' = 'motif'))


# determine shared species and similarity
pr.jaccard <-
  pr.pairs |>
  filter(motif.x < motif.y) |>
  inner_join(pr.species, c('motif.x' = 'motif', 'region'),
             relationship = 'many-to-many') |>
  inner_join(pr.species, c('motif.y' = 'motif', 'region', 'tax.bio')) |>
  count(motif.x, motif.y, name = 'shared') |>
  # add 0 for pairs without overlaps
  right_join(
    pr.pairs |> filter(motif.x < motif.y),
    c('motif.x', 'motif.y')
  ) |>
  mutate_at('shared', replace_na, 0L) |>
  # jaccard similarity on species level
  mutate(jaccard = shared / (n.x + n.y - shared))


p2.jaccard <-
  pr.jaccard |>
  ggplot(aes(jaccard)) +
  stat_ecdf() +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('Jaccard similarity of species shared\nbetween CRSs detected within\nthe same search region') +
  ylab('Empirical cumulative density') +
  theme_bw(18)

################################################################################
################################################################################
# Part 2: Sequence overlaps
################################################################################
################################################################################

# coordinates as ranges
pos.range <-
  pos |>
  select(
    seqnames, start, end, strand,
    motif, tax.bio
  ) |>
  as_granges() |>
  mutate(pos = paste(motif, seqnames, start, end, strand))

regions.range <- regions |> as_granges()

################################################################################
# limit sequence hits to search regions

mask <-
  pos.range |>
  mutate(
    dest = motif |> str_remove('.fna.*$'),
    hit = width
  ) |>
  join_overlap_intersect(regions.range) |>
  filter(dest == region) |>
  filter(width == hit)

pos.range <-
  pos.range |>
  filter(pos %in% mask$pos)
# reduces to 17421 from 17444

################################################################################
# check for overlaps

no.pos <-
  pos.range |>
  as_tibble() |>
  count(motif, name = 'postotal')

overlaps <-
  pos.range |>
  mutate(
    w = width, s = strand
  ) %>%
  join_overlap_intersect(., .) |>
  as_tibble() |>
  # ignore overlap with itself
  filter(pos.x != pos.y) |>
  transmute(
    motif.x, pos.x, len.x = w.x,
    motif.y, pos.y, len.y = w.y,
    overlap = width,
    antisense = s.x != s.y
  ) |>
  left_join(no.pos, c('motif.x' = 'motif')) |>
  left_join(no.pos, c('motif.y' = 'motif'))


# some general informative numbers on the overlap

list(
  no.motifs =  no.pos |> nrow(),
  with.any.overlap = overlaps |>
    select(motif.x, motif.y) |>
    unlist() |>
    unique() |>
    length(),
  with.any.sense.overlap = overlaps |>
    filter(!antisense) |>
    select(motif.x, motif.y) |>
    unlist() |>
    unique() |>
    length(),
  with.self.overlap = overlaps |>
    filter(motif.x == motif.y) |>
    select(motif.x, motif.y) |>
    unlist() |>
    unique() |>
    length(),
  with.self.sense.overlap = overlaps |>
    filter(motif.x == motif.y) |>
    filter(!antisense) |>
    select(motif.x, motif.y) |>
    unlist() |>
    unique() |>
    length()
)
# $no.motifs
# [1] 423
# 
# $with.any.overlap
# [1] 164
# 
# $with.any.sense.overlap
# [1] 153
# 
# $with.self.overlap
# [1] 5
# 
# $with.self.sense.overlap
# [1] 3

################################################################################
# Statistics on the rel overlap lengths

rel.over <-
  overlaps |>
  # focus on the sense overlaps
  filter(! antisense) |>
  # overlap relative to shorter length
  mutate(
    len = pmin(len.x, len.y),
    x = overlap / len
  )

p3.relover <-
  rel.over |>
  # ignore summetric A-B B-A comparisons for general overview of density
  filter(pos.x < pos.y) |>
  ggplot(aes(x)) +
  stat_ecdf() +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_vline(xintercept = .9, color = 'red') +
  geom_hline(yintercept = .15, color = 'red') +
  xlab('Sequence overlap relative to\n shorter sequence length') +
  ylab('Empirical cumulative density') +
  theme_bw(18)

# proportions of alignment with which motif.x sequences overlap above cutoff
aln.prop <-
  rel.over |>
  filter(x >= .9) |>
  # keep symmetric A-B B-A cases, but ignore self hits
  filter(motif.x != motif.y) |>
  count(motif.x, motif.y, postotal.x, postotal.y) |>
  mutate(prop = n / postotal.x)


p4.alnprop <-
  aln.prop |>
  ggplot(aes(prop)) +
  stat_ecdf()  +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_vline(xintercept = .9, color = 'blue') +
  geom_hline(yintercept = .6, color = 'blue') +
  xlab('Proportion of alignment overlapped\n(with rel. overlap ≥ 0.9)') +
  ylab('Empirical cumulative density') +
  theme_bw(18)


################################################################################
# list of motif groups that are potentially redundant to each other

pot.pairs <-
  aln.prop |>
  filter(prop >= .9) |>
  select(motif.x, motif.y) %>%
  # include bi-directional case
  bind_rows(
    .,
    . |>
      mutate(
        tmp = motif.x,
        motif.x = motif.y,
        motif.y = tmp
      ) |>
      select(- tmp)
  ) |>
  unique()

# detect components in 'graph'
# for now a loop approach, because tidygraph is not part of the renv container

# all nodes
nodes <- pot.pairs |> unlist() |> unique()
# length(nodes)
# 85

component.helper <- function(nodes) {
  # starting node
  v0 <- first(nodes)
  # list of nodes it can reach
  reached <- c(v0)
  while(TRUE) {
    # nodes reachable from already reached nodes
    new.reached <-
      union(
        reached,
        pot.pairs |> filter(motif.x %in% reached) |> pull(motif.y)
      )
    # stop if no new nodes have been reached
    if (length(reached) == length(new.reached)) {
      break
    }
    reached <- new.reached
  }
  # ok, no new nodes reached -> one component
  # continue on other nodes, if there are any and combine results
  rest <- setdiff(nodes, reached)
  if(length(rest) > 0) {
    return(c(list(reached), component.helper(rest)))
  } else {
    return(list(reached))
  }
}

components <- component.helper(nodes)

redundant.candidates <-
  components %>%
  set_names(1:length(.)) %>%
  map2(names(.), ~ tibble(g = .y, motif = .x, n = length(.x))) |>
  bind_rows() |>
  arrange(desc(n), g) |>
  # group name sorted by component size
  mutate(
    g2 = lag(g, default = first(g)),
    toggle = g != g2,
    g3 = cumsum(toggle) + 1
  ) |>
  transmute(
    group = paste0('group', g3),
    no.motifs = n,
    motif
  )


################################################################################
# some plot for those candidates

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
    'text',
    3.5, 30,
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
  xlab('CRS in candidate redundant group') +
  ylab('No. candidate redundant group') +
  theme_bw(18)


################################################################################
################################################################################
# Part 3: Check consensus structure for sequence overlapping motifs
################################################################################
################################################################################

# load consensus structure from the stockholm files

consensus.helper <- function(motif) {
  #'data/F_cmfinder/D_search-seqs/K02221_downstream/K02221_downstream.fna.motif.h1_1' |>
  'data/F_cmfinder/D_search-seqs/%s/%s' |>
    sprintf(
      motif |> str_remove('.fna.motif.*$'),
      motif
    ) |>
    read_lines() |>
    # find consensus
    keep(str_detect, '^#=GC SS_cons') |>
    # remove marker and note
    str_remove('^#=GC SS_cons *') |>
    # combine line
    str_c(collapse = '')
}

rc.seq <-
  redundant.candidates |>
  mutate(consensus = motif |>
           map(consensus.helper) |>
           unlist()) |>
  # tidy up dot-bracket
  mutate_at('consensus', str_replace_all, '<', '(') |>
  mutate_at('consensus', str_replace_all, '>', ')') |>
  mutate_at('consensus', str_replace_all, '[^()]', '.')
  
# save for later inspection
rc.seq |>
  write_tsv('foo.tsv')

################################################################################

test.pairs <-
  rc.seq %>%
  left_join(., ., c('group', 'no.motifs')) |>
  filter(motif.x < motif.y)

test.pairs |>
  select(consensus.x, consensus.y) |>
  # remove leading/trailing un-paired bases
  mutate_at(c('consensus.x', 'consensus.y'), str_remove, '^\\.+') |>
  mutate_at(c('consensus.x', 'consensus.y'), str_remove, '\\.+$') |>
  apply(1, str_c, collapse = '\n') |>
  str_c(collapse = '\n') |>
  write_lines('foo.txt')


#singularity exec -B $PWD:$PWD viennarna:2.6.3--py39pl5321h4e691d4_0  bash -c "cat foo.txt | RNAdistance -B > foo2.txt"

rna.dist <-
  'foo2.txt' |>
  read_lines() |>
  keep(str_detect, '^f:') |>
  str_remove('f: ') |>
  as.integer()

################################################################################

test.pairs |>
  mutate(f = rna.dist) |>
  arrange(f) |>
  # write_tsv('stefan.tsv')
  View()

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

cowplot::plot_grid(
  p1.regionbar,
  p2.jaccard,
  p3.relover,
  p4.alnprop,
  p5.components,
  labels = 'AUTO',
  label_size = 20
)
