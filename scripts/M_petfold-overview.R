library(tidyverse)
library(ggplot2)
library(cowplot)
library(furrr)


cpus <- as.integer(unlist(snakemake@threads))
# cpus <- 20
plan(multisession, workers = cpus)

################################################################################
# Parse PETfold output

parser <- function(path) {
  # path <- 'data/M_PETfold/candidates/K00008_upstream.fna.motif.h1_2/output.txt'
  xs <- read_lines(path)
  
  tibble(
    motif = path |>
      dirname() |>
      basename(),
    type = path |>
      dirname() |>
      dirname() |>
      basename(),
    structure = xs |>
      keep(str_starts, pattern = '^PETfold RNA structure:') |>
      str_remove('^PETfold RNA structure:[[:space:]]*'),
    score = xs |>
      keep(str_starts, pattern = '^Score') |>
      str_extract('[0-9.]+$') |>
      as.numeric(),
    ensemble = xs |>
      keep(str_starts, pattern = '^Length-normalized ensemble diversity') |>
      str_extract('[0-9.]+$') |>
      as.numeric()
  )
}


dat <-
  'data/M_PETfold/*/*/output.txt' |>
  Sys.glob() |>
  future_map(parser) |> 
  bind_rows()

################################################################################
# (from J_recall.R)
# Ahead of recall, query Rfam family types for more informative plots

'https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/bacteria.csv' %>%
  read_csv() %>%
  select(Family, Name = `Rfam ID`, type = `RNA type`) -> rfam.types

# keep.type <- c('Cis-reg', 'antisense', 'CRISPR', 'ribozyme', 'rRNA', 'tRNA', 'sRNA')
# 
# rfam.types %>%
#   mutate(
#     type2 = type %>%
#       map(function(i) {
#         keep.type %>%
#           map(~ str_detect(i, .x)) %>%
#           unlist %>%
#           which %>%
#           first -> i
#         keep.type[i]
#       }) %>%
#       unlist
#   ) %>%
#   mutate_at('type2', replace_na, 'other') %>%
#   rename(type.full = type, type = type2) -> rfam.types

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
# (Important also for later, now sto.paths/sto.ss have same order as dat)
sto.paths <- sto.paths[dat$motif]

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
  d <- bit::symdiff(xs, ys)
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
# Compare structures


# extract bps positions
bps.pos <- list(
  petfold = future_map(dat$structure, dotbracket),
  cmfinder = future_map(sto.ss, dotbracket)
)


bp.dist <- with(bps.pos, future_map2(petfold, cmfinder, bp.dist)) |> unlist()

################################################################################
# Combine in single tidy table

petfold <-
  dat |>
  left_join(
    rfam.types |>
      select(motif = Family, class = type),
    'motif'
  ) |>
  mutate_at('class', replace_na, 'CMfinder candidate') |>
  select(- type) |>
  select(class, motif, PETfold.score = score, PETfold.ensemble = ensemble,
         PETfold.structure = structure) |>
  mutate(
    PETfold.bps = bps.pos$petfold |>
      map(length) |>
      unlist(),
    alignment.structure = sto.ss |> unlist(),
    alignment.bps = bps.pos$cmfinder |>
      map(length) |>
      unlist(),
    bp.distance = bp.dist,
    bp.shared = (alignment.bps + PETfold.bps - bp.distance) / 2,
    jaccard = bp.shared / (alignment.bps + PETfold.bps - bp.shared)
  )

write_tsv(petfold, 'data/M_PETfold.tsv')

################################################################################
# Boxplot Ensemble

petfold |>
  left_join(
    petfold |>
      count(class),
    'class'
  ) |>
  mutate(
    xl = ifelse(class == 'CMfinder candidate', 'Predicted novel motif (this study)', ' Bacterial Rfam family') |>
      as.factor(),
    class = sprintf('%s (n=%g)', class, n) |>
      fct_reorder(PETfold.ensemble)
  ) |>
  ggplot(aes(class, PETfold.ensemble, color = xl)) +
  geom_boxplot() +
  scale_y_log10() +
  xlab(NULL) +
  ylab('Length-normalized ensemble diversity') +
  coord_flip() +
  ggsci::scale_color_jama(name = NULL) +
  theme_bw(18) +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical') -> p1


################################################################################
# Boxplot scores

petfold |>
  left_join(
    petfold |>
      count(class),
    'class'
  ) |>
  mutate(
    xl = ifelse(class == 'CMfinder candidate', 'Predicted novel motif (this study)', ' Bacterial Rfam family') |>
      as.factor(),
    class = sprintf('%s (n=%g)', class, n) |>
      fct_reorder(PETfold.score) |>
      fct_rev()
  ) |>
  ggplot(aes(class, PETfold.score, color = xl)) +
  geom_boxplot() +
  # scale_y_log10() +
  xlab(NULL) +
  ylab('Length-normalized PETfold score') +
  coord_flip() +
  ggsci::scale_color_jama(name = NULL) +
  theme_bw(18) +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical') -> p2


################################################################################
################################################################################
# helper for loading reliability matrix
rel.mat <- function(path) {
  # path <- 'data/M_PETfold/candidates/K00008_upstream.fna.motif.h1_2/reliabilities.txt'
  lines <-
    path |>
    read_lines()
  
  # Expected alignment size
  i <- lines[[1]] |> as.integer()
  
  # Load matrix
  lines[-1] |>
    I() |>
    read_delim(delim = ' ', col_names = FALSE, skip_empty_rows = TRUE) |>
    as.matrix() -> res
  colnames(res) <- NULL
  
  # Last column in NA because all lines with a space -> remove
  assertthat::assert_that(res[, ncol(res)] |> is.na() |> all())
  res <- res[, - ncol(res)] 
  
  # remove last row
  res <- res[-nrow(res), ]
  assertthat::are_equal(dim(res), c(i, i))
  return(res)
}


# Load reliability matrices
'data/M_PETfold/*/*/reliabilities.txt' |>
  Sys.glob() %>%
  set_names(. |> dirname() |> basename()) |>
  future_map(safely(rel.mat)) -> mats

# Exclude those that did not have enough sequences for PETfold to run
# mats |>
#   map('error') |>
#   discard(is.null) |>
#   length()
# 146
mats |>
  map('result') |>
  discard(is.null) -> mats
  

mats |>
  map(isSymmetric) |>
  unlist() |>
  table()
# TRUE 
# 1350 

################################################################################
################################################################################
# Ratio for scores on bp vs total sum

bps <-
  petfold |>
  with(set_names(PETfold.structure, motif)) |>
  map(dotbracket) |>
  map(~ invoke(rbind, .x))


petfold |>
  pull(motif) |>
  map(~ mats[[.x]][bps[[.x]]] |> sum()) |>
  unlist() -> struct.sum

petfold |>
  pull(motif) |>
  map(function(i) {
    x <- mats[[i]]
    sum(x[upper.tri(x)])
  }) |>
  unlist() -> upper.sum

petfold |>
  transmute(
    class, motif,
    PETfold.score, PETfold.ensemble,
    struct.sum = struct.sum,
    upper.sum = upper.sum,
    ratio = struct.sum / upper.sum
  ) |>
  left_join(
    petfold |>
      count(class),
    'class'
  ) |>
  mutate(
    xl = ifelse(class == 'CMfinder candidate', 'Predicted novel motif (this study)', ' Bacterial Rfam family') |>
      as.factor(),
    class = sprintf('%s (n=%g)', class, n) |>
      fct_reorder(ratio) |>
      fct_rev()
  ) |>
  ggplot(aes(class, ratio, color = xl)) +
  geom_boxplot() +
  xlab(NULL) +
  ylab('Base-pairing reliability ratio\nPaired positions vs total sum') +
  coord_flip() +
  ggsci::scale_color_jama(name = NULL) +
  theme_bw(18) +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical') -> p3


################################################################################
################################################################################
plot_grid(p2, p1, p3, nrow = 1,
          labels = 'AUTO', label_size = 20)

ggsave('data/M_PETfold-stats.png', width = 28, height = 10, dpi = 400)

################################################################################
