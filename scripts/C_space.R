library(tidyverse)
library(parallel)
library(ape)

path.trees <- 'data/C_phylo/*/bootstrap-consensus.tree'

# make sure script output is placed in log file
# log <- file(unlist(snakemake@log), open="wt")
# sink(log)

# the numer of cores to use
# cpus <- as.integer(unlist(snakemake@threads))
cpus <- 16

cl <- makeForkCluster(cpus)

################################################################################

xs <- path.trees %>%
  Sys.glob() %>%
  set_names(. %>% dirname %>% basename) %>%
  map(read.tree)

################################################################################

# Remove multi-gene entries, only keep species info on tip
despec <- function(t) {
  # t <- xs$K00135
  tibble(
    gene = t$tip.label,
    spec = str_remove(gene, '\\.[^.]*$'),
    med = apply(cophenetic.phylo(t), 1, median)
  ) %>%
    group_by(spec) %>%
    slice_min(med, n = 1, with_ties = FALSE) %>%
    ungroup() -> dat
  t2 <- keep.tip(t, dat$gene)
  t2$tip.label <- dat$spec
  return(t2)
}

# small test example in which A.1.b should be removed
assertthat::are_equal(
  read.tree(text = '((A.1.a:1,A.1.b:1.5):2,(B.2.a:1,C.3.a:1):3);') %>%
    despec() %>%
    unroot %>%
    reorder.phylo("postorder"),
  read.tree(text = '(A.1:3,(B.2:1,C.3:1):3);') %>%
    unroot %>%
    reorder.phylo("postorder")
)

# despec input trees
xs2 <- xs %>% map(despec)

################################################################################
# For a pair, focus on shared species and compute distance
dedist <- function(x) {
  # a <- xs2$K00135 ; b <- xs2$K02699
  a <- x[[1]] ; b <- x[[2]]
  a <- xs2[[a]] ; b <- xs2[[b]]
  shared <- intersect(a$tip.label, b$tip.label)
  
  a2 <- keep.tip(a, shared)
  b2 <- keep.tip(b, shared)
  
  #  Kuhner and Felsenstein scoring as described in
  # Kuhner M. K. & Felsenstein J. 1994. Simulation comparison of phylogeny
  # algorithms under equal and unequal evolutionary rates. Molecular Biology
  # and Evolution 11: 459–468.
  list(
    shared = length(shared),
    dist = dist.topo(a2, b2, 'score')
  )
}

# All pairwise comparisons
crossing(a = names(xs), b = names(xs)) -> tasks

# res <- parRapply(cl, head(tasks), safely(dedist))
res <- parRapply(cl, tasks, safely(dedist))

res %>%
  map(function(x) {
    if (is.null(x$result)) {
      tibble(
        shared = NA_real_,
        dist = NA_real_
      )
    } else {
      tibble(
        shared = as.numeric(x$result$shared),
        dist = as.numeric(x$result$dist)
      )
    }
  }) %>%
  bind_rows() -> res2
    
# bind_cols(head(tasks), res2)
pds <- bind_cols(tasks, res2)

################################################################################

# Inspect to what extend the observed distances depend on the number of shared
# species

cat('Correlation test of no. shared species to topology distance')
with(pds, cor.test(shared, dist)) %>% print()

pds %>%
  ggplot(aes(shared, dist)) +
  geom_hex(bins=100) +
  scale_fill_viridis_c() +
  geom_smooth(method = 'lm', se = FALSE, color = 'red') +
  xlab('No. shared species') +
  ylab('Pair-wise tree topology distance') +
  theme_bw(18)
  

################################################################################
# Explore the space of trees accordint to a PCoA of the distance matrix

# Issue: some values are NA, impute with average
pds %>%
  drop_na() %>%
  group_by(a) %>%
  summarize(avg = mean(dist)) %>%
  ungroup %>%
  right_join(pds, 'a') %>% 
  mutate(d2 = ifelse(is.na(dist), avg, dist)) -> pds.imputed

# Buld matrix
pds.imputed %>%
  select(a, b, d2) %>%
  spread(b, d2) -> foo
foo %>%
  select(- a) %>%
  as.matrix() %>%
  magrittr::set_rownames(foo$a) -> pds.mat
# assure symmetry, due to average numerically not perfectly symmetric in individual
# cases, but ignore for now
# assertthat::assert_that(isSymmetric(pds.mat))


# pds.mat[upper.tri(pds.mat)] %>%
#   hist

# helper to plot MDS scaled distances and highligh potential outliers
mds.helper <- function(mat, outlier = 3) {
  pds.scl <- mat %>%
    as.dist() %>%
    cmdscale()
  
  
  pds.scl %>%
    magrittr::set_colnames(c('x', 'y')) %>%
    as_tibble(rownames = 'OG') %>%
    mutate(
      dx = x - mean(x),
      dy = y - mean(y),
      sdx = abs(dx) / sd(x),
      sdy = abs(dy) / sd(y),
      mx = sdx > outlier,
      my = sdy > outlier,
      m = mx | my,
      l = ifelse(m, OG, NA_character_)
    ) -> dat
  dat %>%
    ggplot(aes(x, y, label = l)) +
    geom_point(size = 3, alpha = 0.5) +
    ggrepel::geom_label_repel(size = 7) +
    xlab('MDS 1') +
    ylab('MDS 2') +
    theme_bw(18) -> p

  list(dat, ggExtra::ggMarginal(p, type = 'hist'))
}

foo <- mds.helper(pds.mat)
pds.scl <- foo[[1]]
pds.scl.mds <- foo[[2]]

pds.scl.mds

################################################################################
# Query for mot extreme outlier for assessment

pds.scl %>%
  mutate(msd = pmax(sdx, sdy)) %>%
  slice_max(msd, n = 1, with_ties = FALSE) %>%
  pull(OG) -> extremeout

xs[[extremeout]] %>%
  ladderize() %>%
  plot(show.tip.label = FALSE)

################################################################################
################################################################################
# Follow the idea of
# Mai, Uyen, and Siavash Mirarab.
# “TreeShrink: Fast and Accurate Detection of Outlier Long Branches in
#   Collections of Phylogenetic Trees.”
# BMC Genomics 19, no. S5 (May 2018): 272.
# https://doi.org/10.1186/s12864-018-4620-2.

# Out of curiosity keep track of the length of splits
# -> might provide an additional level of support

# For comparison compute overall split distance distribution
split.dists <- function(x) {
  # x <- xs2$K05592
  x <- reorder.phylo(unroot(x), "postorder")
  nx <- length(x$tip.label)
  # in the reordered unrooted tree, all internal nodes have indices
  # >= nx + 2
  # splits are
  bps <- ape:::bipartition2(x$edge, nx)
  # and the boot split distances are those corresponding to entries in which
  # the edge table's 2nd column has an internal node
  mask <- between(x$edge[, 2], nx + 2, nx + length(bps))
  tmp <- tibble(
    splitdist =  x$edge.length[ mask ],
    # the 'left split' sizes, first entry is full tree
    size = map(bps[ -1 ], length) %>% unlist
  ) %>%
    mutate(
      # check the 'right split' and take minimum
      size = pmin(size, length(x$tip.label) - size)
    ) 
  # Also include distances of the tips, though not a split, individual
  # species might also not belong to the tree
  # those are edges starting at tips, which are the lowest index
  tibble(
    splitdist =  x$edge.length[ x$edge[, 2] <= length(x$tip.label) ],
    size = 1
  ) %>%
    bind_rows(tmp)
}

parLapply(cl, xs, split.dists) %>%
  map2(names(.), ~ mutate(.x, OG = .y)) %>%
  bind_rows() -> overall.splits


# Compute tree diameter inflation of splits
# Here, the algorithm loop structure is similar to the
# split distances computed by
# ape:::.dist.topo.score

# tree diamter
helper.diameter <- function(x) {
  dist.nodes(x) %>% max
}

# For a set of tip indices, remove corresponding tips
# and compute diamter
helper.trim.diam <- function(x, ix, drop = TRUE) {
  if (drop) {
    x <-  drop.tip(x, ix)
  } else {
    x <-  keep.tip(x, ix)
  }
  helper.diameter(x)
}

dia.ratio <- function(x) {
  x.dia <- helper.diameter(x)
  
  # Start loop similar to split.dists
  x <- reorder.phylo(unroot(x), "postorder")
  nx <- length(x$tip.label)
  # in the reordered unrooted tree, all internal nodes have indices
  # >= nx + 2
  # splits are
  bps <- ape:::bipartition2(x$edge, nx)
  # and the boot split distances are those corresponding to entries in which
  # the edge table's 2nd column has an internal node
  
  # the 'left split' sizes, first entry is full tree, exclude
  bps <- bps[ -1 ]
  
  # diameters in left split
  bps %>%
    map(helper.trim.diam, x = x) %>%
    unlist -> l.dia
  # diameters in right split
  bps %>%
    map(helper.trim.diam, x = x, drop = FALSE) %>%
    unlist -> r.dia
  
  # For splits, check impact on diamter
  mask <- between(x$edge[, 2], nx + 2, nx + length(bps) + 1)
  tibble(
    splitdist =  x$edge.length[ mask ],
    # the 'left split' sizes, first entry is full tree
    l.size = map(bps, length) %>% unlist,
    l.dia = l.dia,
    r.size = nx - l.size,
    r.dia = r.dia,
    bps = bps
  )  -> tmp
    
  # determine which side of the split is the largest part remaining
  tmp %>%
    mutate(
      left = l.size > r.size,
      trim.size = ifelse(left, l.size, r.size),
      trim.dia = ifelse(left, l.dia, r.dia)
    ) %>%
    transmute(
      bps, trim.size, trim.dia, splitdist,
      # the bps indicate what was removed, so make sure to remember
      # that if the result is from the right, then those should be kept
      keep.tip = ! left
    ) -> tmp2
    
  # Also include distances of the tips, though not a split, individual
  # species might also not belong to the tree
  # those are edges starting at tips, which are the lowest index
  tibble(
    splitdist =  x$edge.length[ x$edge[, 2] <= length(x$tip.label) ],
    trim.size = nx - 1,
    keep.tip = FALSE,
    bps = as.list(1:nx),
    trim.dia = map(1:nx, helper.trim.diam, x = x) %>% unlist
  ) %>%
    bind_rows(tmp2) -> tmp3
  
  # determine which split maximize the diamter shrink factor
 tmp3  %>%
   mutate(shrink = x.dia / trim.dia) %>%
   slice_max(shrink) %>%
   # in case of multiple maxia, take the one with most tips remaining
   slice_max(shrink) %>%
   # and if then still, the next best one with max splitdist
   slice_max(splitdist, with_ties = FALSE)
}

res <- parLapply(cl, xs, safely(dia.ratio))
parLapply(cl, xs, dia.ratio) %>%
  map2(names(.), ~ mutate(.x, OG = .y)) %>%
  bind_rows() -> shrinks


################################################################################
################################################################################
################################################################################
################################################################################



stopCluster(cl)