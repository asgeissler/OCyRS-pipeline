# Load sequences of motif alignments that were used to scan the genome for homologs
# Determine via string lookup location in genomes

library(tidyverse)
library(furrr)

# in.fdr <- 'data/I_fdr.tsv'
# config.cutoff <- 10

in.fdr <- unlist(snakemake@input)
config.cutoff <- unlist(snakemake@config[['fdr_candidates']])

out.seq <- unlist(snakemake@output[['seq']])
out.pos <- unlist(snakemake@output[['pos']])

cpus <- as.integer(unlist(snakemake@threads))
# cpus <- 20
plan(multisession, workers = cpus)


###############################################################################
# Load sequences of alignments

worker <- function(path, motif) {
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
    mutate_at('seq', str_remove_all, '[^A-Z]') %>%
    # keep track of position of lines
    mutate(r = 1:n()) %>%
    # combine for genes with info in multiple lines
    # in order of apperance in file
    arrange(gene, r) %>%
    group_by(gene) %>%
    summarize(stockholm.seq = str_c(seq, collapse = '')) %>%
    ungroup %>%
    mutate(tax.bio = str_remove(gene, '\\.[^.]*$')) %>%
    select(tax.bio, stockholm.seq) %>%
    mutate(motif = motif)
}

in.fdr %>%
  read_tsv() %>%
  filter(RNAphylo.fdr <= config.cutoff) %>%
  mutate(path = sprintf(
    'data/F_cmfinder/D_search-seqs/%s/%s',
    region, motif
  )) %>%
  select(path, filename = motif) %>%
  with(future_map2(path, filename, worker)) %>% 
  bind_rows() -> motif.seqs
  

###############################################################################
# Find locations of all aligned sequences following the
# Efficient genome searching tutorial, section
# "Finding an arbitrary nucleotide pattern in an entire genome"
# https://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/GenomeSearching.pdf
# Idea: reverse complement pattern before search

motif.seqs$tax.bio %>%
  unique %>%
  future_map(function(i) {
    sprintf('data/A_representatives/%s/genome.fna.gz', i) %>%
      Biostrings::readDNAStringSet() -> genome
    
    motif.seqs %>%
      filter(tax.bio == i) %>%
      mutate(row = 1:n()) %>%
      mutate(strand = '+') -> fwd
    
    fwd %>%
      mutate(
        stockholm.seq = stockholm.seq %>%
          Biostrings::DNAStringSet() %>%
          Biostrings::reverseComplement() %>%
          as.character(),
        strand = '-'
      ) %>%
      bind_rows(fwd) -> patterns
    
    
    patterns %>%
      group_by_all() %>%
      do(search = {
        grp <- .
        Biostrings::vmatchPattern(grp$stockholm.seq, genome) %>%
          as.data.frame() %>%
          transmute(start, end, chr = names(genome)[group])
      }) %>%
      ungroup %>%
      unnest(search) %>%
      select(motif, chr, start, end, strand)
  }) %>%
  bind_rows() -> motif.seqs.pos


###############################################################################

write_tsv(motif.seqs, out.seq)
write_tsv(motif.seqs.pos, out.pos)
