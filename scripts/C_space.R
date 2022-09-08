library(tidyverse)
library(furrr)
library(ape)

path.raw.trees <- 'data/C_phylo/*/bootstrap-consensus.tree'
path.shrunk.trees <- 'data/C_shrink/*/output.tree'
path.ref.trees <- 'reference-trees/*.tree'

# make sure script output is placed in log file
# log <- file(unlist(snakemake@log), open="wt")
# sink(log)

# the numer of cores to use
# cpus <- as.integer(unlist(snakemake@threads))
cpus <- availableCores()
plan(multisession, workers = cpus) 

################################################################################
# glob lookup all tree files
raw <- path.raw.trees %>%
  Sys.glob() %>%
  set_names(. %>% dirname %>% basename)

shrunk <- path.shrunk.trees %>%
  Sys.glob() %>%
  set_names(. %>% dirname %>% basename)

ref <- path.ref.trees %>%
  Sys.glob() %>%
  set_names(. %>% basename %>% fs::path_ext_remove())

################################################################################
# load trees
  
raw.trees    <- future_map(c(raw, ref), read.tree)
shrunk.trees <- future_map(c(shrunk, ref), read.tree)

################################################################################

# De-duplicate species entries.
# The tip lapels are taxid.bioproject.gene
# or taxid.gcf.name / taxid.name for the reference trees
# Per species (taxid) keep gene that is on average closest to alll other
# leaves in the tree (mean cophenetic distance).
# Afterward, only keep taxid in labels
despec <- function(t) {
  tibble(
    gene = t$tip.label,
    taxid = str_remove(gene, '\\..*$'),
    avgdist = apply(cophenetic.phylo(t), 1, mean)
  ) %>%
    group_by(taxid) %>%
    # only keep the one of lowest avgdist,
    # but exactly one, if with ties choose arbriatrily
    slice_min(avgdist, n = 1, with_ties = FALSE) %>%
    ungroup() -> dat
  t2 <- keep.tip(t, dat$gene)
  t2$tip.label <- dat$taxid
  return(t2)
}

# small test example
assertthat::are_equal(
  read.tree(text = '((A.1.a:1,A.1.b:1.5):2,(B.2.a:1,C.3.a:1):3);') %>%
    despec() %>%
    unroot %>%
    reorder.phylo("postorder"),
  read.tree(text = '(A:3,(B:1,C:1):3);') %>%
    unroot %>%
    reorder.phylo("postorder")
)

# despec input trees
raw.despec.trees    <- future_map(raw.trees, despec)
shrunk.despec.trees <- future_map(shrunk.trees, despec)

################################################################################
# For a pair, of species de-duplicated trees, compute the topology distance
# on the set of shared species
dedist <- function(trees, a, b) {
  a <- trees[[a]] ; b <- trees[[b]]
  # focus on shared nodes
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
assertthat::are_equal(names(raw.despec.trees),
                      names(shrunk.despec.trees))
tasks <- crossing(a = names(raw.despec.trees),
                  b = names(shrunk.despec.trees)) %>%
  filter(a < b)

# suppress erroneous warnings of
# "unexpectedly generated random numbers without specifying argument 'seed'."
# seems to be a bug in furrr
f.raw <- safely(partial(dedist, trees = raw.despec.trees))
f.shrunk <- safely(partial(dedist, trees = shrunk.despec.trees))

suppressWarnings({
  res.raw <- future_pmap(tasks, f.raw)
  res.shrunk <- future_pmap(tasks, f.shrunk)
})

################################################################################
# Collect distances as table

mk.tbl <- function(i) {
  # i <- res.raw
  i %>%
    map('error') %>%
    map(is.null) %>%
    unlist -> ok.mask
  table(ok.mask)
  
  i[ok.mask] %>%
    map('result') %>%
    future_map(function(j) {
      j$dist <- as.numeric(j$dist)
      as_tibble(j)
    }) %>%
    bind_rows() -> dat
  
  dat2 <- bind_cols(tasks[ok.mask, ], dat)
}

res.raw.tbl <- mk.tbl(res.raw)
res.shrunk.tbl <- mk.tbl(res.shrunk)

################################################################################
# Inspect to what extend the observed distances depend on the number of shared
# species

cat('Correlation test of no. shared species to topology distance')
cat('Raw trees')
with(res.raw.tbl, cor.test(shared, dist)) %>% print()
cat('Shrunk trees')
with(res.shrunk.tbl, cor.test(shared, dist)) %>% print()

################################################################################
# Build matrices
# note: The distance is possible to fail if not enough species were shared
# idea: substitute by max observed distance

mk.mat <- function(i) {
  # i <- res.shrunk.tbl
  # Create empty matrix
  ns <- i %>%
    select(a, b) %>%
    unlist %>%
    unique
  i.mat <- matrix(NA_real_, nrow = length(ns), ncol = length(ns)) %>%
    magrittr::set_rownames(ns) %>%
    magrittr::set_colnames(ns)
  # Fill values both triangles
  i.mat[cbind(i$a, i$b)] <- i$dist
  i.mat[cbind(i$b, i$a)] <- i$dist
  # the distance per gene at which outlier start
  dm <- apply(i.mat, 1, function(x) {
    quantile(x, .75, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE) 
  })
  
  dm.vec <-rep(dm, each = length(ns))
  dm.vec2 <-rep(dm, times = length(ns))
  # The vector can be used to make a matrix by:
  # matrix(dm.vec, nrow = length(ns), ncol = length(ns))
  # thus the vectors correspond to the transposition of each other
  # => use to compute local max observed outlier distance between a pair
  map2(dm.vec, dm.vec2, c) %>%
    map(max) %>%
    unlist -> dm.dat
  
  dm.mat <-matrix(dm.dat, nrow = length(ns), ncol = length(ns))
  assertthat::assert_that(isSymmetric(dm.mat))
  
  mask <- is.na(i.mat)
  i.mat[mask] <- dm.dat[mask]
  
  # Fill diagonal, should be after the averaging to not influence impute value
  diag(i.mat) <- 0
  
  # assure symmetry
  assertthat::assert_that(isSymmetric(i.mat))
  return(i.mat)
}

raw.mat <- mk.mat(res.raw.tbl)
shrunk.mat <- mk.mat(res.shrunk.tbl)

################################################################################
# Explore the space of trees according to a PCoA of the distance matrix

list(
  'Full trees' = raw.mat,
  'Outlier pruned' = shrunk.mat
) %>%
  map(as.dist) %>%
  map(cmdscale) %>%
  map(magrittr::set_colnames, c('MDS1', 'MDS2'))  %>%
  map(as_tibble, rownames = 'tree') %>%
  map2(names(.), ~ mutate(.x, mode = .y)) %>%
  bind_rows() -> mds

mds %>%
 muate(ref = !str_detect(tree, '^K[0-9]*$')) -> dat
dat %>% filter(ref) -> mds.ref
dat %>%
  ggplot(aes(MDS1, MDS2, color = ref)) +
  ggsci::scale_color_jama() +
  geom_point(size = 3, alpha = 0.5) +
  geom_point(size = 3, alpha = 0.5, data = mds.ref) +
  ggrepel::geom_label_repel(
    aes(label = tree), data = mds.ref,
    size = 7,
    # nudge_x = 10, nudge_y = -5,
    force = 500,
    alpha = 0.7, show.legend = FALSE) +
  facet_wrap(~ mode, scales = 'free') +
  theme_bw(18)  +
  theme(legend.position = 'hide')

################################################################################
