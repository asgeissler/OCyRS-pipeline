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

all.genes.ranges <-
  all.genes |>
  select(seqnames = tax.bio.chr, start, end, strand, gene = tax.bio.gene) |>
  as_granges()

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
    a85 = xs.ranges$K05585@seqnames@values,
    a06 = xs.ranges$K02906@seqnames@values
  )
)
# - 178 pairs on same chromosome
# - Genes without pairs (either not in genome or different chromosome)
#   K05585 - 16
#   K02906 - 14

################################################################################
#Q3 distance between anchors

gene.dist <- plyranges::add_nearest_distance(xs.ranges$K05585, xs.ranges$K02906) %>%
  filter(!is.na(distance))

gene.dist$distance %>% summary
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     245     392     707  265100    1174 5336622  

################################################################################
# Follow-up question:
# What is the distance and motif present distribution per species?


tax <-
  'data/A_representatives/taxonomy.tsv' |>
  read_tsv()


list( 
  xs.ranges %>%
    # extract tax.bio for genes
    map(~ GenomicRanges::seqnames(.x)) %>%
    map(as.character) %>%
    map(str_remove, '\\.[^.]*$') %>%
    map2(names(.), ~ tibble(
      tax.bio = .x,
      x = sprintf('Gene-copies of %s', .y),
      n = 1
    )) %>%
    bind_rows() %>%
    group_by(tax.bio, x) %>%
    summarize(n = sum(n)) %>%
    ungroup %>%
    spread(x, n),
  
  tibble(
    tax.bio = GenomicRanges::seqnames(gene.dist) |>
      as.character() |>
      str_remove('\\.[^.]*$'),
    'Distance gene pair' = gene.dist$distance
  ) %>%
    group_by(tax.bio) %>%
    summarize_all(min) %>%
    ungroup,
  
  aln.pos.range %>%
    as_tibble() %>%
    mutate(
      tax.bio = seqnames %>%
        as.character() %>%
        str_remove('\\.[^.]*$'),
      n = 1
    ) %>%
    group_by(tax.bio, motif) %>%
    summarize(n = sum(n)) %>%
    ungroup %>%
    mutate(x = paste('No. motif pos', motif)) %>%
    select(- motif) %>%
    spread(x, n)
) %>%
  purrr::reduce(.f = full_join, by = 'tax.bio') %>%
  left_join(tax, 'tax.bio') -> followup
  
View(followup)

# - 195 genome have at least one copy of either ortholog gene
# - 1 has 2 gene copies of the K02906 and a different genome has 2 copies for K05585
# - 178 genomes have both genes on the same contig/choromosome
# 
# - 4 genomes have only a single copy of either but not both (xor) gene
# - 2 of these 4 genomes have the K05585_upstream.fna.motif.h1_3
# 
# - the remaining 13 genome have both genes but on different contains -> unknown distance
# - one of the 13 has K02906_upstream.fna.motif.h1_2


followup %>% 
  drop_na(`Distance gene pair`) %>%
  mutate(
    foo = case_when(
    !is.na(`No. motif pos K02906_upstream.fna.motif.h1_2`) & !is.na(`No. motif pos K05585_upstream.fna.motif.h1_3`) ~ 'both motifs',
    !is.na(`No. motif pos K02906_upstream.fna.motif.h1_2`) &  is.na(`No. motif pos K05585_upstream.fna.motif.h1_3`) ~ 'either motif',
     is.na(`No. motif pos K02906_upstream.fna.motif.h1_2`) & !is.na(`No. motif pos K05585_upstream.fna.motif.h1_3`) ~ 'either motif',
    TRUE ~ 'neither motif'
    )
  ) %>%
  # count(foo)
#  both motifs      48
#  either motif      4
#  neither motif   126
  ggplot(aes(foo, `Distance gene pair`)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, color = 'blue' , alpha = .7) +
  scale_y_log10() +
  xlab(NULL) +
  theme_bw(18)


################################################################################
# Follow-up question:
# Those with both: Is there a correlation between adjacency and distance?

with.both <- as.character(seqnames(gene.dist))


xs.ranges %>%
  set_names(NULL) %>%
  invoke(.f = plyranges::pair_nearest) %>%
  as_tibble %>%
  transmute(
    seqnames = granges.x.seqnames,
    start = pmin(granges.x.end, granges.y.end) + 1,
    end = pmax(granges.x.start, granges.y.start) - 1,
    strand = '*'
  ) %>%
  as_granges -> min.gaps



join_overlap_intersect(min.gaps, all.genes.ranges) %>%
  seqnames() %>%
  as.character() %>%
  str_remove('\\.[^.]*$') %>%
  unique -> not.adj

followup %>%
  select(
    tax.bio,
    `Distance gene pair`,
  ) %>%
  drop_na %>%
  mutate(x = ifelse(
    tax.bio %in% not.adj,
    "Not adjacent",
    "Adjacent"
  )) %>%
  ggplot(aes(x, `Distance gene pair`)) +
  geom_boxplot() +
  scale_y_log10()
