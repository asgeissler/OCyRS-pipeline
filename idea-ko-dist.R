# There were two potentially reverse-complementary motifs
# ? Are the search regions maybe overlapping (same motif)
# ? Orientation of the search anchors
# ? Distances between them

library(tidyverse)
library(plyranges)

# The motifs of interest were
xs <- c(
  'K05585_upstream.fna.motif.h1_3',
  'K02906_upstream.fna.motif.h1_2'
)
################################################################################

# Q1 overlap between motif's alignment sequences
aln.pos <-
  'data/K_motif-tax-pos.tsv' |>
  read_tsv() |>
  filter(motif %in% xs)

aln.pos.range <-
  aln.pos |>
  select(seqnames, start, end, strand, motif) |>
  as_granges() |>
  mutate(s = strand, w = width)

n1 <- filter(aln.pos.range, motif == xs[[1]]) |> length()
n2 <- filter(aln.pos.range, motif == xs[[2]]) |> length()

join_overlap_intersect(
  filter(aln.pos.range, motif == xs[[1]]),
  filter(aln.pos.range, motif == xs[[2]])
) |>
  mutate(
    shared = width,
    jaccard = shared / (w.x + w.y - shared),
    orient = ifelse(s.x == s.y, 'sense', 'anti-sense')
  )  -> foo

foo$orient %>% unique
# only anti-sense
foo$jaccard %>% ecdf %>% plot

n3 <- foo |>
  filter(jaccard >= .9) |>
  length()

n3 / (n1 + n2 - n3)
# 78% of positions very similar (>= 90%) with anti-sense overlaps

################################################################################
#Q2 orientation of the search anchors

xs.anchors <-
  xs |>
  str_remove('_upstream.fna.motif.*$')

# Gene ids used in search region (phylogeny filtering etc)
xs.genes <-
  sprintf('data/C_shrink/%s/output.tree', xs.anchors) |>
  map(ape::read.tree) |>
  map(~. $tip.label) |>
  set_names(xs.anchors)

# all gene coordinates
all.genes <- read_tsv('data/A_representatives/genes.tsv.gz')

xs.ranges <-
  xs.genes |>
  map(~ tibble(tax.bio.gene = .x)) |>
  map(~ semi_join(all.genes, .x, 'tax.bio.gene')) |>
  map(select, seqnames = tax.bio.chr, start, end, strand) |>
  map(as_granges)

xs.ranges %>% map(length)
# - Genes used for analysis (passed MSA+phylogeny filtering):
#   K05585 - 195
#   K02906 - 193

venn::venn(
  list(
    a5 = xs.ranges$K05585@seqnames@values,
    a6 = xs.ranges$K02906@seqnames@values
  )
)
# - 178 pairs on same chromosome
# - Genes without pairs (either not in genome or different chromosome)
#   K05585 - 16
#   K02906 - 14

################################################################################
#Q3 distance between anchors

foo <- plyranges::add_nearest_distance(xs.ranges$K05585, xs.ranges$K02906) 

foo$distance %>% summary
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#     245     392     707  265100    1174 5336622      16 