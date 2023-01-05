# Purpose: Collect cmstat output into single table

library(tidyverse)


c(
  'data/I_cmstat/*/*.txt',
  'data/I_cmstat-rfam/*.txt'
) %>%
  map(Sys.glob) %>%
  unlist %>%
  map(read_lines) %>%
  map(discard, str_detect, pattern = '^#') %>%
  unlist -> xs

cols <- c(
  'idx', 'name', 'accession',
  'nseq', 'eff_nseq',
  'clen', 'W',
  'bps', 'bifs',
  'model',
  'entropy.cm', 'entropy.hmm'
)

xs %>%
  # remove leading spaces
  str_remove('^ *') %>%
  # split files on whitespace
  str_split(pattern = ' +') %>%
  # set colnames
  map(set_names, cols) %>%
  # collect as table
  bind_rows() %>%
  # keep Rfam idea instead of name for matching
  mutate(
    name = ifelse(accession == '-', name, accession)
  ) %>%
  # nicer names
  select(
    'motif' = name,
    nseq, eff_nseq,
    consensus_residues_len = clen,
    expected_max_hit_len = W,
    bps,
    bifurcations = bifs,
    rel_entropy_cm = entropy.cm,
    rel_entropy_hmm = entropy.hmm
  ) %>%
  mutate_at(vars(- motif), as.numeric) -> xs


write_tsv(xs, 'data/I_cmstat.tsv')