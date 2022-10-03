# Purpose:
# Extract from SISSIz text output the shuffled alignment
# Export alignment FastA
# Double check that seqid/dinucleotide content/GC are (more ore less) the same

library(tidyverse)
library(furrr)
library(ape)
library(Biostrings)

in.full <- 'data/D_search-seqs-aln/*.fna.gz'
in.filtered <- 'data/E_search-filtered/*.fna.gz'
in.shuffled <- 'data/E_search-shuffled/*.txt'

out <- 'data/E_search-filtered-stat.tsv'
# side-effect: Convert SISSIz txt to fna.gz

# make sure script output is placed in log file
log <- file(unlist(snakemake@log), open="wt")
sink(log)

# cpus <- 16
cpus <- as.integer(unlist(snakemake@threads))
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
# Parse SISSIz

helper <- function(i) {
  # Read put skip the first long line with the tree
  dat <- read_lines(i, skip = 1)
  # write a temporary clustal file with a valid header
  tmp <- tempfile()
  header <- 'CLUSTAL X (1.81) multiple sequence alignment'
  dat <- c(header, dat)
  write_lines(dat, tmp)
  # read the temporary file
  dat <- readDNAMultipleAlignment(tmp, format = 'clustal')
  # cleanup temp file
  file.remove(tmp)
  return(dat)
}

in.shuffled %>%
  Sys.glob() %>%
  set_names(. %>% basename %>% str_remove('.txt$')) %>%
  future_map(safely(helper)) -> alns.shuffled

# Collect those runs for which SISSIz succeeded
alns.shuffled %>%
  map('error') %>%
  map(is.null) %>%
  unlist() -> mask
table(mask)
alns.shuffled[mask] %>%
  map('result') -> alns.shuffled


# order of rows in alignments should be the same
# But: SISSIz truncated the names
# Manually rename to full name: Needed for subsequent CMfinder screen and
# scoring of random models!
alns.filtered[names(alns.shuffled)] %>%
  map(rownames) -> names.shouldbe

map2(
  alns.shuffled, names.shouldbe,
  function(x, y) {
    rownames(x) <- y
    x
  }
) -> alns.shuffled

# Make a StringSet without the gaps
alns.shuffled %>%
  future_map(function(i) {
  # map(function(i) {
    i <- as(i, 'DNAStringSet')
    i[] <- str_remove_all(i, '-')
    i
  }) -> alns.shuffled.nogaps
      

################################################################################
# Write out alignments in fasta form (needed for cmfinder)

alns.shuffled %>%
  future_map2(names(.), ~ writeXStringSet(
    as(.x, 'DNAStringSet'),
    sprintf('%s/%s.aln.fna.gz', 'data/E_search-shuffled', .y),
    compress = TRUE
  ))

alns.shuffled.nogaps %>%
  future_map2(names(.), ~ writeXStringSet(
    .x,
    sprintf('%s/%s.fna.gz', 'data/E_search-shuffled', .y),
    compress = TRUE
  ))

################################################################################
# General statistics

helper <- function(i, j) {
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
  future_map2(names(.), helper) %>%
  bind_rows() -> stat.full

alns.filtered %>%
  future_map2(names(.), helper) %>%
  bind_rows() -> stat.filtered

alns.shuffled %>%
  future_map2(names(.), helper) %>%
  bind_rows() -> stat.shuffled

stat.full %>%
  full_join(
    stat.filtered,
    suffix = c('', '.filtered'),
    'region'
  ) %>%
  full_join(
    stat.shuffled,
    suffix = c('', '.shuffled'),
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
  alns.filtered[todo$region],
  alns.shuffled[todo$region],
  di.helper
) %>%
  map2(names(.), ~tibble(region = .y, dinucleotide.freq.cor = .x)) %>%
  bind_rows() -> di.res

################################################################################
# combine stats and export

stat.res %>%
  left_join(di.res, 'region') -> res

write_tsv(res, out)

# Report final correlation of gc and seqid
print('Correlation of overall GC content before/after filtering')
print(with(res, cor.test(gc, gc.filtered)))
print('Correlation of overall GC content filtered vs shuffled')
print(with(res, cor.test(gc.shuffled, gc.filtered)))
print('Correlation of overall avg. seq. id content before/after filtering')
print(with(res, cor.test(avg.seqid, avg.seqid.filtered)))
print('Correlation of overall avg. seq. id content filtered vs shuffled')
print(with(res, cor.test(avg.seqid.shuffled, avg.seqid.filtered)))
print(with(res, plot(avg.seqid.shuffled, avg.seqid.filtered)))

print('Average di-nucleotide freq. correlation')
res$dinucleotide.freq.cor %>% summary %>% print
print('Standard deviation')
res$dinucleotide.freq.cor %>% sd


print('DONE.')
