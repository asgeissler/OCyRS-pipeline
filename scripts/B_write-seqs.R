# Purpose: For each candidate OG, write out the peptide sequneces

library(tidyverse)
library(Biostrings)

path.tbl <- unlist(snakemake@input[['tbl']])
path.seqs <- unlist(snakemake@input[['seqs']])

path.out <- unlist(snakemake@output)

# path.tbl <- 'data/B_OGs.tsv'
# path.seqs <- 'data/A_progenomes/proteins.faa.gz'

# path.out <- 'data/B_OGs'

###############################################################################
# The candidate OGs

tbl <- read_tsv(path.tbl)

###############################################################################
# Per genome the peptide sequences

file.path(path.seqs, '*/proteins.fasta.gz') %>%
  Sys.glob() %>%
  set_names(. %>% dirname %>% basename) %>%
  map(readAAStringSet) -> ps

###############################################################################
# Per KO write out corersponding sequences

dir.create(path.out)
for(ko in unique(tbl$term)) {
  # ko <- 'K00003'
  print(ko)
  tbl %>%
    filter(term == ko) %>%
    select(tax.bio, tax.bio.gene) %>%
    group_by(tax.bio) %>%
    do(j = ps[[dplyr::first(.$tax.bio)]][.$tax.bio.gene]) %>%
    ungroup() %>%
    pull(j) %>%
    purrr::reduce(c) -> xs
  
  writeXStringSet(xs, sprintf('%s/%s.faa.gz', path.out, ko),
                  compress = TRUE)
}


###############################################################################
###############################################################################