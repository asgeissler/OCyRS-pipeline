# Purpose: Collect cmstat output into single table

library(tidyverse)


'data/I_cmstat/*/*.txt' %>%
  Sys.glob() %>%
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
  str_remove('^ *') %>%
  str_split(pattern = ' +') %>%
  map(set_names, cols) %>%
  bind_rows() %>%
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