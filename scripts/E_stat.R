# Purpose:
# Check for gappy regions in alignments. These should be removed before
# the SISSIz shuffeling

library(tidyverse)
library(furrr)
library(ape)
library(Biostrings)

in.full <- 'data/D_search-seqs-aln/*.fna.gz'
in.filtered <- 'data/E_search-filtered/*.fna.gz'

out <- 'data/E_search-filtered-stat.tsv'

# make sure script output is placed in log file
log <- file(unlist(snakemake@log), open="wt")
sink(log)

cpus <- 16
# cpus <- as.integer(unlist(snakemake@threads))
plan(multisession, workers = cpus) 

################################################################################
# Load alignments

alns.full <- in.full %>%
  Sys.glob() %>%
  set_names(. %>% basename %>% str_remove('.fna.gz$')) %>%
  future_map(readDNAMultipleAlignment)

alns.filtered <- in.filtered %>%
  Sys.glob() %>%
  set_names(. %>% basename %>% str_remove('.fna.gz$')) %>%
  future_map(readDNAMultipleAlignment)

################################################################################
# General statistics

helper <- function(alns, i, j) {
  # i <- alns[[1]]
  x <- as.DNAbin(i)
  
  #Part1: Average Sequence identity
  # Pairwise distances
  sid <- dist.gene(x, pairwise.deletion = TRUE, method = 'percentage')
  sids <- sid[upper.tri(sid)]
  # dist to similarity
  sids <- (1 - sids) * 100
  sid <- mean(sids, na.rm = TRUE)
  
  #Part2: GC content
  gc <- GC.content(x)
  
  tibble(
    region = j,
    avg.seqid = sid,
    gc = gc
  )
}

alns.full %>%
  future_map2(names(.), partial(helper, alns = .)) %>%
  bind_rows() -> stat.full

alns.filtered %>%
  future_map2(names(.), partial(helper, alns = .)) %>%
  bind_rows() -> stat.filtered

full_join(
  stat.full,
  stat.filtered,
  suffix = c('.full', '.filtered'),
  'region'
) -> stat.res

################################################################################
# Compare di-nucleotide contents

di.helper <- function(x, y) {
  # i <- 'K00003_downstream'
  # x <- alns.full[[i]]
  # y <- alns.filtered[[i]]
  
  x %>%
    as('DNAStringSet') %>%
    dinucleotideFrequency() %>%
    unlist -> xs
  y %>%
    as('DNAStringSet') %>%
    dinucleotideFrequency() %>%
    unlist -> ys
  
  cor.test(xs, ys)$estimate
}

# only work on those in which a filtered region exist
stat.res %>%
  drop_na %>%
  select(region) -> todo

future_map2(
  alns.full[todo$region],
  alns.filtered[todo$region],
  di.helper
) %>%
  map2(names(.), ~tibble(region = .y, dinucleotide.freq.cor = .x)) %>%
  bind_rows() -> di.res

################################################################################
# combine stats and export

stat.res %>%
  left_join(di.res, 'region') -> res

write_tsv(seq, out)

# Report final correlation of gc and seqid
print('Correlation of overall GC content before/after filtering')
print(with(res, cor.test(gc.full, gc.filtered)))
print('Correlation of overall avg. seq. id content before/after filtering')
print(with(res, cor.test(avg.seqid.full, avg.seqid.filtered)))
