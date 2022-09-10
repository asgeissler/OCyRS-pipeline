library(tidyverse)
library(furrr)
library(ape)

# path.raw.trees <- 'data/C_phylo/*/bootstrap-consensus.tree'
# path.shrunk.trees <- 'data/C_shrink/*/output.tree'
# path.ref.trees <- 'reference-trees/*.tree'

path.raw.trees <- unlist(snakemake@input[['raw']])
path.shrunk.trees <- unlist(snakemake@input[['shrunk']])
path.ref.trees <- unlist(snakemake@input[['refs']])

out.raw <- unlist(snakemake@output[['raw']])
out.shrunk <- unlist(snakemake@output[['shrunk']])
out.mds <- unlist(snakemake@output[['mds']])
out.mds.fig <- unlist(snakemake@output[['mdsfig']])


# make sure script output is placed in log file
log <- file(unlist(snakemake@log), open="wt")
sink(log)

# the numer of cores to use
# cpus <- availableCores()
cpus <- as.integer(unlist(snakemake@threads))
plan(multisession, workers = cpus) 

################################################################################
# (glob lookup all tree files during development only)
raw <- path.raw.trees %>%
  # Sys.glob() %>%
  set_names(. %>% dirname %>% basename)

shrunk <- path.shrunk.trees %>%
  # Sys.glob() %>%
  set_names(. %>% dirname %>% basename)

ref <- path.ref.trees %>%
  # Sys.glob() %>%
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
  a2 <- keep.tip(a, shared) %>% unroot()
  b2 <- keep.tip(b, shared) %>% unroot()
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

write_tsv(res.raw.tbl, out.raw)
write_tsv(res.shrunk.tbl, out.shrunk)

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
  # i <- res.raw.tbl
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
  diag(i.mat) <- 0
  i.mat[cbind(i$a, i$b)] <- i$dist
  i.mat[cbind(i$b, i$a)] <- i$dist
  # Replace na 
  # i.mat[is.na(i.mat)] <- quantile(i$dist, .75, na.rm = TRUE) +
  #   1.5 * IQR(i$dist, na.rm = TRUE)
  # i.mat[is.na(i.mat)] <- max(i$dist)
  i.mat[is.na(i.mat)] <- mean(i$dist)
  
  # Fill diagonal, should be after the averaging to not influence impute value
  
  # assure symmetry
  assertthat::assert_that(isSymmetric(i.mat))
  return(i.mat)
}

# raw.mat <- mk.mat(res.raw.tbl)
shrunk.mat <- mk.mat(res.shrunk.tbl)
# scale by #shared
res.shrunk.tbl %>%
  mutate(dist = dist * shared) %>%
  mk.mat() -> shrunk2.mat

################################################################################
# Explore the space of trees according to a PCoA of the distance matrix

list(
  'Full trees' = raw.mat,
  'Outlier pruned' = shrunk.mat,
  'Scaled by no. shared species' = shrunk2.mat
) %>%
  map(as.dist) %>%
  map(cmdscale) %>%
  map(magrittr::set_colnames, c('MDS1', 'MDS2'))  %>%
  map(as_tibble, rownames = 'tree') %>%
  map2(names(.), ~ mutate(.x, mode = .y)) %>%
  bind_rows() -> mds

mds %>%
 mutate(ref = !str_detect(tree, '^K[0-9]*$')) -> dat
dat %>% filter(ref) -> mds.ref
dat %>%
  ggplot(aes(MDS1, MDS2, color = ref)) +
  ggsci::scale_color_jama(name = 'Reference tree') +
  geom_point(size = 3, alpha = 0.5) +
  geom_point(size = 3, alpha = 0.5, data = mds.ref) +
  ggrepel::geom_label_repel(
    aes(label = tree), data = mds.ref,
    size = 7,
    # nudge_x = 10, nudge_y = -5,
    # force = 500,
    force = 100,
    alpha = 0.7, show.legend = FALSE) +
  facet_wrap(~ mode, scales = 'free') +
  xlim(-5, 10) +
  ylim(-6, 6) +
  theme_bw(18)  +
  theme(legend.position = 'bottom')

write_tsv(dat, out.mds)
ggsave(out.mds.fig, dpi = 400, width = 16, height = 9)
################################################################################
