library(tidyverse)
library(ggplot2)
library(cowplot)
library(furrr)

# Assume exists, the motifs
#  data/J_novel/potentially-novel-motifs.tsv
# Rfam files, eg with implicit dependency via
# data/G_rfam-cmsearch.tsv.gz

cpus <- as.integer(unlist(snakemake@threads))
# cpus <- 20
plan(multisession, workers = cpus)

################################################################################
# Parse RNAalifold output

parser <- function(path) {
  # path <- 'data/M_RNAalifold/candidates/K02906_upstream.fna.motif.h2_5'
  # Alignment length
  aln.len <-
    path |>
    file.path('alifold.out') |>
    read_lines(n_max = 1) |>
    str_remove('^.*length of alignment ') |>
    as.integer()
  
  # MFE + ensemble
  xs <-
    path |>
    file.path('stdout.txt') %>%
    read_lines()
  
  xs2 <-
    xs[[length(xs) - 1]] |>
    str_trim() |>
    str_remove(' =.*') |>
    str_split(' \\{') |> 
    unlist()
    
  mfe.struct <- xs2[[1]]  
  mfe.val <- xs2[[2]]
  
  xs2 <-
    xs[[length(xs)]] |>
    str_trim() |>
    str_remove('frequency of mfe structure in ensemble ') |>
    str_split('; ensemble diversity') |>
    unlist()
  
  mfe.freq <- xs2[[1]]
  ensemble <- xs2[[2]]
  
  # Parse base-pairing probabilities
  bps <-
    path |>
    file.path('alidot.ps') |>
    read_lines() |>
    keep(str_detect, '[ul]box$') |>
    keep(str_detect, '^[0-9]') |>
    str_remove('^.*hsb ')
  bps <-
    tibble(
      foo = str_remove(bps, ' [ul]box$'),
      upper = str_detect(bps, ' ubox$')
    ) |>
    separate(foo, c('i', 'j', 'sqrt.prob'), sep = ' ') |>
    mutate(across(c(i, j), as.integer)) |>
    mutate(across(sqrt.prob, as.numeric)) |>
    mutate(
      tmp = i,
      i = ifelse(upper, i, j),
      j = ifelse(upper, j, tmp)
    )
  
  bps.mat <- matrix(0, nrow = aln.len, ncol = aln.len)
  bps.mat[with(bps, cbind(i, j))] <- bps$sqrt.prob
  
  # results
  list(
    alignment.length = aln.len,
    mfe.structure = mfe.struct,
    mfe.kcal.mol = mfe.val,
    mfe.freq = mfe.freq,
    ensemble.diversity = ensemble,
    bp.probs = bps.mat
  )
}

################################################################################

paths <-
  'data/M_RNAalifold/*/*' |>
  Sys.glob()

dat <-
  paths |>
  future_map(parser)

# Tidy up results as table

ali.tbl <-
  dat |>
  map(`[`, c('alignment.length', 'mfe.structure', 'mfe.kcal.mol',
        'mfe.freq', 'ensemble.diversity')) |>
  bind_rows() |>
  mutate_all(str_trim) |>
  mutate(across( - 'mfe.structure', as.numeric)) |>
  mutate(
    model = paths |> basename(),
    type = paths |> dirname() |> basename()
  )

bp.probs <-
  dat |>
  map('bp.probs') |>
  set_names(basename(paths))

################################################################################
# load CMfinder predicted consensus structure

sto.paths <-
  c(
    'data/F_cmfinder/D_search-seqs/*/*.fna.motif.*',
    'data/G_rfam-bacteria-seeds/*.sto'
  ) |>
  map(Sys.glob) |>
  unlist() %>%
  set_names(. |> basename() |> str_remove('.sto$'))

# Only keep needed ones
sto.paths <- sto.paths[ali.tbl$model] 

# Extract predicted structure
sto.ss <-
  sto.paths |>
  future_map(function(path) {
    # path <- 'data/G_rfam-bacteria-seeds/RF02470.sto'
    path |>
      read_lines() |>
      keep(str_detect, '^#=GC SS_cons') |>
      str_remove('^#=GC SS_cons *') |>
      str_c(collapse = '')
  })

# Simplify notation (see infernal userguide)
# http://eddylab.org/infernal/Userguide.pdf

sto.ss <-
  sto.ss |>
  map(str_replace_all, '<', '(') |>
  map(str_replace_all, '>', ')') |>
  map(str_replace_all, '[^()]', '.')

################################################################################
# Helper for extracting basepairs and distance

dotbracket <- function(xs) {
  res <- list()
  todo <- c()
  for (i in seq(1, nchar(xs))) {
    x <- str_sub(xs, i, i)
    if (x == '(') {
      todo <- c(todo, i)
    } else if (x == ')') {
      j <- todo[[length(todo)]]
      todo <- todo[- length(todo)]
      res <- c(res, list(c(j, i)))
    }
  }
  res
}

assertthat::are_equal(dotbracket(''), list())
assertthat::are_equal(dotbracket('(.)'), list(c(1, 3)))
assertthat::are_equal(dotbracket('((.)(.))'), list(c(2, 4), c(5, 7), c(1, 8)))

my.sort <- function(xs) {
  if (length(xs) == 0) {
    xs
  } else {
    xs2 <- map(xs, sort)
    xs2.first <- map(xs2, 1) |> unlist()
    xs2[order(xs2.first)]
  }
}

assertthat::are_equal(my.sort(list()), list())
assertthat::are_equal(
  my.sort(list(
    c(2, 3),
    c(5, 4),
    c(1, 6)
  )),
  list(
    c(1, 6),
    c(2, 3),
    c(4, 5)
  )
)

bp.dist <- function(xs, ys) {
  xs <- my.sort(xs)
  ys <- my.sort(ys)
  d <- symdiff(xs, ys)
  length(d)
}

assertthat::are_equal(
  bp.dist(dotbracket('((..))'),
          dotbracket('((..))')),
  0
)
assertthat::are_equal(
  bp.dist(dotbracket('()..()'),
          dotbracket('((..))')),
  4
)


################################################################################
# Compare structures of CMfinder and RNAalifold

cmp.tbl <-
  ali.tbl |>
  mutate(cmfinder.structure = sto.ss[model] |> unlist())

# extract bps positions
bps.pos <- with(cmp.tbl, list(
  mfe = future_map(mfe.structure, dotbracket),
  cmfinder = future_map(cmfinder.structure, dotbracket)
))

cmp.tbl <-
  cmp.tbl |>
  mutate(bp.dist = with(bps.pos, future_map2(mfe, cmfinder, bp.dist)) |> unlist())


################################################################################
# First plots

# Distinguish Rfam with/without Cyanobacteria hits
hit.cyano <-
  'data/G_rfam-cmsearch.tsv.gz' |>
  read_tsv() |>
  select(rf) |>
  unique()
cmp.tbl2 <-
  cmp.tbl |>
  semi_join(hit.cyano, c('model' = 'rf')) |>
  mutate(type = 'Rfam, with CMsearch hit in Cyanobacteria') |>
  bind_rows(cmp.tbl)

# Ensemble diversity
p1 <-
  cmp.tbl2 |>
  ggplot(aes(ensemble.diversity / (alignment.length ** 2), color = type)) +
  stat_ecdf(size = 1.2) +
  ggsci::scale_color_jama(name = NULL) +
  ylab('Empirical cumulative density') +
  xlab('Ensemble diversity / Alignment length ^ 2') +
  theme_bw(16) +
  theme(legend.position = 'bottom')

# MFE
p2 <-
  cmp.tbl2 |>
  ggplot(aes(mfe.kcal.mol / (alignment.length ** 2), color = type)) +
  stat_ecdf(size = 1.2) +
  ggsci::scale_color_jama(name = NULL) +
  ylab('Empirical cumulative density') +
  xlab('Free energy [kcal / mol] / Alignment length ^ 2') +
  theme_bw(16) +
  theme(legend.position = 'bottom')

# MFE Freq
p3 <-
  cmp.tbl2 |>
  ggplot(aes(mfe.freq, color = type)) +
  stat_ecdf(size = 1.2) +
  ggsci::scale_color_jama(name = NULL) +
  ylab('Empirical cumulative density') +
  xlab('MFE structure frequency') +
  theme_bw(16) +
  theme(legend.position = 'bottom')

# CMfinder vs MFE
p4 <-
  cmp.tbl2 |>
  ggplot(aes(bp.dist / (alignment.length ** 2), color = type)) +
  stat_ecdf(size = 1.2) +
  ggsci::scale_color_jama(name = NULL) +
  ylab('Empirical cumulative density') +
  xlab('CMfinder/Rfam annotated structure vs MFE\nbp distance / Alignment length ^ 2') +
  theme_bw(16) +
  theme(legend.position = 'bottom')

plot_grid(p1, p2, p3, p4,
          labels = 'AUTO', label_size = 16)

ggsave('foo.png', width = 16, height = 14)

################################################################################

cmp.tbl2 |>
  mutate(i = 1:n()) |>
  filter(model == 'K02906_upstream.fna.motif.h2_5') |>
  pull(i) -> i
mat <- bp.probs[[i]]
f <- upper.tri
f2 <- lower.tri

mat[f2(mat)] <- mat[f(mat)]

hist(1 - apply(1-mat, 1, prod))


################################################################################

