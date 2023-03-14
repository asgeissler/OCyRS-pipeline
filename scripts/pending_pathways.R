left_join(
  potential.novel %>%
    count(category, name = 'potential.novel.motifs'),
  potential.novel %>%
    left_join(motifs) %>%
    filter(str_detect(chr, '^32049\\.')) %>%
    select(motif, category) %>%
    unique %>%
    count(category, name = 'motifs.in.PCC7002')
)


potential.novel %>%
  left_join(motifs) %>%
  filter(str_detect(chr, '^32049\\.')) %>%
  select(motif, category) %>%
  head