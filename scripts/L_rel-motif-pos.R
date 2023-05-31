# Match RNA motif sequence alignment positions relative to search regions
# The script structure writes a file similar to "J_motif-aln-seq-pos.tsv" 
# But: The positions are relative to the search region, NOT the genome
# Then: Compare rel. positions to CMsearch hits

library(tidyverse)
library(furrr)

cpus <- as.integer(unlist(snakemake@threads))
# cpus <- 30
plan(multisession, workers = cpus)

################################################################################
# Query which motifs and search regions were of interest

interest.regions <-
  'data/L_regions-todo.txt' |>
  read_lines() %>%
  tibble(region = .)

interest.motifs <-
  'data/J_novel/potentially-novel-motifs.tsv' |>
  read_tsv()

################################################################################
# Alignment sequences

all.seqs <-
  'data/J_motif-aln-seqs.tsv.gz' |>
  read_tsv() |>
  semi_join(interest.motifs, 'motif') |>
  mutate(region = str_remove(motif, '\\.fna\\.motif.*$'))

################################################################################
# Adapted from 'J_aln-seq-pos.R' relative position lookup

interest.regions |>
  pull(region) |>
  future_map(function(i) {
    # i <- 'K00008_upstream'
    sprintf('data/D_search-seqs/%s.fna.gz', i) %>%
      Biostrings::readDNAStringSet() -> region
    
    all.seqs |>
      filter(region == i) |> 
      select(stockholm.seq) |>
      unique() |>
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
        Biostrings::vmatchPattern(grp$stockholm.seq, region) %>%
          as.data.frame() |>
          transmute(start, end, chr = names(region)[group])
      }) %>%
      ungroup %>%
      unnest(search) %>%
      select(stockholm.seq, chr, start, end, strand) |>
      mutate(region = i)
  }) %>%
  bind_rows() -> motif.rel.aln.pos

all.seqs |>
  left_join(
    motif.rel.aln.pos |>
      mutate(tax.bio = str_remove(chr, '\\.[^.]*$')),
    c('stockholm.seq', 'region', 'tax.bio')
  ) |>
  select(- stockholm.seq, - tax.bio) -> rel.aln.pos

rel.aln.pos |>
  write_tsv('data/L_motif-aln-rel-pos.tsv.gz')

################################################################################
# Load CMsearch results
# Script in part adatped from 'G_combine_search.R'



# names of columns in file
ns <- c('target name', 'accession', 'query name', 'accession2',
        'mdl', 'mdl from', 'mdl to', 'seq from', 'seq to', 'strand',
        'trunc', 'pass', 'gc', 'bias', 'score',
        'E-value', 'inc', 'description of target')

parse <- function(path) {
  # path <- 'data/L_cmsearch-regions/K00008_upstream/K00008_upstream.fna.motif.h1_2.txt'
  ls <- read_lines(path)
  
  # spaces in second line are field breaks
  pos <- str_locate_all(ls[2], ' ')[[1]]
  assertthat::are_equal(pos[, 'start'], pos[, 'end'])
  
  # Extract per non-comment line the fields
  ls %>%
    discard(str_detect, '^#') %>%
    map(function(i) {
      map2(
        c(1, pos[, 'start']),
        c(pos[, 'start'], nchar(ls[2])),
        ~ str_sub(i, .x, .y)
      ) %>%
        map(str_trim) %>%
        unlist()
    }) %>%
    invoke(.f = rbind) %>%
    as_tibble() %>%
    # set column names as from first line comment
    magrittr::set_colnames(ns) %>%
    # remember folder name is tax_bio / genome id
    mutate(
      region = path %>%
        dirname %>%
        basename
    )
}

# Parse all files
# 'data/L_cmsearch-regions/*/*.txt' %>% # Too many files for both '*' :(
'data/L_cmsearch-regions/*/' |>
  Sys.glob() |>
  map(paste0, '*.txt') |>
  map(Sys.glob) |>
  unlist() |>
  future_map(parse) %>%
  invoke(.f = bind_rows) -> tbl

# more concise names
select(
  tbl,
  region, chr = `target name`,
  start = `seq from`, end = `seq to`, strand,
  motif = `query name`,
  score, evalue = `E-value`, gc
) -> cmsearch.rel.pos

cmsearch.rel.pos |>
  write_tsv('data/L_cmsearch-rel-pos.tsv.gz')

################################################################################
################################################################################