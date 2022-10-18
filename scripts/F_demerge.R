# CMfinder merges motifs yet still reports the merged motifs, for instance
# among the repeated motifs, two are redundant and skew the true number of
# detected motifs:
# ...
# K00228_upstream.fna.motif.h2_4
# K00228_upstream.fna.motif.h1_5
# K00228_upstream.fna.motif.h1_5.h2_4
#
# Purpose: This script removes these merged errenesouly removed motifs
#          from the file list

library(tidyverse)

# path.in <- 'data/predump/motifs.txt'
# path.in <- 'data/F_cmfinder/motifs.txt'

path.in <- unlist(snakemake@input)

out.path <- unlist(snakemake@output)

# make sure script output is placed in log file
log <- file(unlist(snakemake@log), open="wt")
sink(log)

################################################################################
# Load data

# Load raw list
path.in %>%
  map(read_lines) %>%
  unlist() -> raw.list
raw.list <- tibble(path = raw.list) %>%
  mutate(
    filename = basename(path),
    dir = path %>%
      dirname %>%
      basename
  )

raw.list %>%
  separate(filename, c('region', 'helices'),
           sep = '.fna.motif.', remove = FALSE) %>%
  mutate(submotifs = str_count(helices, 'h')) %>%
  separate_rows(helices, sep = '\\.') -> dat

################################################################################
# Detect which submotifs are part of other motif

# Overlap on shared helices
dat %>%
  left_join(dat, c('dir', 'region', 'helices')) %>%
  filter(filename.x != filename.y) %>%
  count(
    dir,
    filename.x, submotifs.x, 
    filename.y, submotifs.y
  ) %>%
  ungroup %>%
  rename(shared = n) %>%
  # files in which all its motifs are contained in another
  filter(submotifs.x == shared) %>%
  select(dir, filename = filename.x) -> merged

################################################################################
# Remove the merged motifs

raw.list %>%
  anti_join(merged, c('dir', 'filename')) -> demerged

################################################################################
# some summarizing tables

print('% of motifs after de-merging')
left_join(
  raw.list %>%
    count(dir) %>%
    rename(raw.nr.of.motifs = n),
  demerged %>%
    count(dir) %>%
    rename(demerged.motifs = n),
  'dir'
) %>%
  mutate(proportion = demerged.motifs / raw.nr.of.motifs * 100) %>%
  print()

################################################################################
# Write out single table of all relevant motifs

demerged %>%
  separate(filename, c('region', 'helices'),
           sep = '.fna.motif.', remove = FALSE) -> res

assertthat::assert_that(
  res %>%
    pull(path) %>%
    map(file.exists) %>%
    unlist %>%
    all
)

write_tsv(res, out.path)