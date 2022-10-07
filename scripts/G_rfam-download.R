library(tidyverse)

# The output dir for the pipeline
out <- unlist(snakemake@output[[1]])
rfv <- unlist(snakemake@params[['rfamv']])

# Query for bacterial Rfams
'https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/bacteria.csv' %>%
  read_csv() %>%
  pull(Family) -> bac

# Load raw cm file
'http://ftp.ebi.ac.uk/pub/databases/Rfam/%s/Rfam.cm.gz' %>%
  sprintf(rfv) %>%
  read_lines() -> dat

# split at '//' into chunks, BUT only every second one due to HMMER3 blocks
newpos <- which(dat == '//')
splitpos <- newpos[seq(2, length(newpos) - 1, 2)]

# intervals for each entry
map2(c(1, splitpos + 1), c(splitpos - 1, length(dat) - 1), seq) %>%
  # partition
  map(function(i) dat[i]) -> entries

# extract Rfam ids
entries %>%
  map(~ .x[[3]]) %>%
  unlist() %>%
  str_remove('^ACC[[:space:]]+') -> ne

assertthat::are_equal(length(ne), length(entries))

# add ending '//' to each entry
entries %>%
  map(~ c(.x, '//')) %>%
  set_names(ne) -> nes

# Write out bacterial cm models
if (!dir.exists(out)) {
  dir.create(out)
}
nes[bac] %>%
  map2(names(.), ~ write_lines(.x, sprintf('%s/%s.cm', out, .y)))

fs::file_touch(sprintf('%s/download.done', out))
