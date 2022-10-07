library(tidyverse)

# The output dir for the pipeline
out <- unlist(snakemake@output[[1]])
rfv <- unlist(snakemake@params[['rfamv']])


# Query for bacterial Rfams
'https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/bacteria.csv' %>%
  read_csv() %>%
  pull(Family) -> bac

# The seeds are used for the scoring steps
# -> exclude the too large ribosomal RNA
rRNA <- c('RF02541', 'RF00177')
bac <- setdiff(bac, rRNA)


# Load raw seed file
'http://ftp.ebi.ac.uk/pub/databases/Rfam/%s/Rfam.seed.gz' %>%
  sprintf(rfv) %>%
  read_lines() -> dat

splitpos <- which(dat == '//')

# intervals for each entry
map2(c(1, splitpos + 1), c(splitpos - 1, length(dat) - 1), seq) %>%
  # partition
  map(function(i) dat[i]) -> entries

# extract Rfam ids
entries %>%
  map(~ .x[[3]]) %>%
  unlist() %>%
  str_remove('#=GF AC[[:space:]]+') -> ne

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
  map2(names(.), ~ write_lines(.x, sprintf('%s/%s.sto', out, .y)))

fs::file_touch(sprintf('%s/download.done', out))
