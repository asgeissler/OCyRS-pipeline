# Adapted from J_aln-seq-pos.R
# Convert Stockholm alignment to FastA (including gaps)

library(tidyverse)
library(furrr)

# Assume exists, the motifs
#  data/J_novel/potentially-novel-motifs.tsv
# Rfam files, eg with implicit dependency via
# data/G_rfam-cmsearch.tsv.gz

cpus <- as.integer(unlist(snakemake@threads))
# cpus <- 20
plan(multisession, workers = cpus)


###############################################################################
# Helper for reading/writing

worker <- function(inpath, outpath) {
  # inpath <- 'data/G_rfam-bacteria-seeds/RF00001.sto'
  # inpath <- 'data/F_cmfinder/D_search-seqs/K00003_upstream/K00003_upstream.fna.motif.h1_1'
  # Load the relative positions from the Stockholm motif file
  inpath %>%
    # read in Stockholm
    read_lines %>%
    # Keep lines with nucleotide info
    keep(str_detect, '^[A-Za-z0-9]') %>%
    tibble(row = .) %>%
    # separate gene id and sequence
    separate(row, c('gene', 'seq'), sep = ' +') %>%
    # !Difference to J_aln-seq-pos.R -> Keep gaps etc
    # mutate_at('seq', str_replace_all, 'U', 'T') %>%
    # mutate_at('seq', str_remove_all, '[^A-Z]') %>%
    # keep track of position of lines
    mutate(r = 1:n()) %>%
    # combine for genes with info in multiple lines
    # in order of appearance in file
    arrange(gene, r) %>%
    group_by(gene) %>%
    summarize(stockholm.seq = str_c(seq, collapse = '')) %>%
    ungroup %>%
    with(set_names(stockholm.seq, gene)) |>
    Biostrings::RNAStringSet() |>
    Biostrings::writeXStringSet(outpath)
}

###############################################################################
# 1. Work on all bacterial Rfam families

'data/G_rfam-bacteria-seeds/*.sto' |>
  Sys.glob() -> xs

dir.create('data/M_alignments-Rfam/')
xs |>
  basename() |>
  str_remove('\\.sto') |>
  sprintf(fmt = 'data/M_alignments-Rfam/%s.fna') -> xs2

future_map2(xs, xs2, worker)

###############################################################################
# 2. Work on the ~400 candidate motifs

'data/J_novel/potentially-novel-motifs.tsv' |>
  read_tsv() |>
  pull(motif) -> mots

sprintf(
  'data/F_cmfinder/D_search-seqs/%s/%s',
  str_remove(mots, '.fna.motif.*'),
  mots
) -> xs

dir.create('data/M_alignments-candidates/')
xs2 <- sprintf('data/M_alignments-candidates/%s.fna', basename(xs))

future_map2(xs, xs2, worker)

           