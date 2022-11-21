# Purpose: Aggregate cmsearch of motifs
library(tidyverse)
library(furrr)

in.cmsearch <- 'data/I_cmsearch/*/*.txt'

cpus <- as.integer(unlist(snakemake@threads))
plan(multisession, workers = cpus) 
################################################################################

# names of columns in file
ns <- c('tax.bio.chr', 'accession', 'motif', 'accession2',
        'mdl', 'mdlfrom', 'mdlto', 'seqfrom', 'seqto', 'strand',
        'trunc', 'pass', 'gc', 'bias', 'score',
        'E-value', 'inc', 'description of target')

parse <- function(path) {
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
      tax_bio = path %>%
        dirname %>%
        basename
    )
}

# Parse all files
'data/I_cmsearch/*/*.txt' %>%
  Sys.glob() %>%
  future_map(parse) %>%
  invoke(.f = bind_rows) -> tbl

################################################################################

# more concise names
select(
  tbl,
  tax_bio,
  tax.bio.chr,
  start = seqfrom, end = seqto, strand,
  name = motif, score, evalue = `E-value`, gc
) -> tbl2

# save 
write_tsv(tbl2, unlist(snakemake@output))

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
