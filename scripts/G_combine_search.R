library(tidyverse)
library(furrr)


cpus <- as.integer(unlist(snakemake@threads))
plan(multicore, workers = cpus) 
################################################################################

# names of columns in file
ns <- c('target name', 'accession', 'query name', 'accession2',
  'mdl', 'mdl from', 'mdl to', 'seq from', 'seq to', 'strand',
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
'data/G_rfam-cmsearch/*/*.txt' %>%
  Sys.glob() %>%
  future_map(parse) %>%
  invoke(.f = bind_rows) -> tbl

################################################################################

# more concise names
select(
  tbl,
  tax_bio, chr = `target name`,
  start = `seq from`, end = `seq to`, strand,
  rf = accession2, name = `query name`,
  score, evalue = `E-value`, gc
) -> tbl2

# save 
write_tsv(tbl2, unlist(snakemake@output))
