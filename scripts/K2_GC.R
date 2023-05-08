# (Parts are adapted from J_aln-seq-pos.R and E_stat.R)
# For the novel motif predictions, compute the GC of the alignments
# (Neither CMfinder/RNAphylo/Rscape explicitly reported it, and so war has
#  only been computed on the search region)

library(tidyverse)
library(furrr)

library(ape)

in.novel <- 'data/J_novel/potentially-novel-motifs.tsv'


# cpus <- as.integer(unlist(snakemake@threads))
cpus <- 20
plan(multisession, workers = cpus)


###############################################################################

novel <- read_tsv(in.novel) %>%
  mutate(
    region = str_remove(motif, '\\.fna.*$'),
    path = sprintf(
      'data/F_cmfinder/D_search-seqs/%s/%s',
      region, motif
    )
  )

###############################################################################
# Load alignments

worker.aln <- function(path) {
  # path <- 'data/F_cmfinder/D_search-seqs/K00003_upstream/K00003_upstream.fna.motif.h1_1'
  # Load the relative positions from the Stockholm motif file
  path %>%
    # read in Stockholm
    read_lines %>%
    # Keep lines with nucleotide info
    keep(str_detect, '^[0-9]') %>%
    tibble(row = .) %>%
    # separate gene id and sequence
    separate(row, c('gene', 'seq'), sep = ' +') %>%
    mutate_at('seq', str_replace_all, 'U', 'T') %>%
      # change relative to J_aln-seq-pos.R to also keep gaps)
    mutate_at('seq', str_remove_all, '[^A-Z-]') %>%
    # keep track of position of lines
    mutate(r = 1:n()) %>%
    # combine for genes with info in multiple lines
    # in order of apperance in file
    arrange(gene, r) %>%
    group_by(gene) %>%
    summarize(stockholm.seq = str_c(seq, collapse = '')) %>%
    ungroup %>%
    mutate(tax.bio = str_remove(gene, '\\.[^.]*$')) %>%
    select(tax.bio, stockholm.seq) -> xs
  
  with(xs, set_names(stockholm.seq, tax.bio)) %>%
    Biostrings::DNAMultipleAlignment()
}

novel %>%
  pull(path) %>%
  future_map(worker.aln) -> alns
  

###############################################################################
# Compute GC content
  
worker.GC <- function(i, j) {
  # i <- alns[[1]]
  x <- as.DNAbin(i)
  GC.content(x)
}

alns %>%
  future_map(worker.GC) %>%
  unlist -> alns.gc

novel %>%
  select(- region, - path) %>%
  mutate(GC = alns.gc) %>%
  write_tsv('data/K2_motifs.tsv')
