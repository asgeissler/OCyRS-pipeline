# Collect scores of R-scape and CMfinder from the various files and
# combine as a single table.

library(tidyverse)
library(furrr)

out.scores <- unlist(snakemake@output)

# implicit dependecies, not passed directly or it would be too much
motifs <- 'data/F_cmfinder/motifs-demerged.tsv'
cm.rfam <- 'data/H_scores/cmscore-Rfam/*.txt'
cm.motifs <- 'data/H_scores/cmscore/*/*.txt'
rscape <- 'data/H_scores/rscape/*/*/stdout.txt'

# cpus <- 8
cpus <- as.integer(unlist(snakemake@threads))
plan(multisession, workers = cpus)

################################################################################
################################################################################
################################################################################
# collect hmmpair and RNAphylo for the Rfam families

cm.rfam %>%
  Sys.glob() %>%
  set_names(. %>% basename %>% fs::path_ext_remove()) %>%
  future_map2(names(.), function(i, j) {
    # i <- 'data/H_scores/cmscore-Rfam/RF00001.txt'
    i <- read_lines(i)
    tibble(
      RNAphylo = keep(i, str_starts, 'RNAphylo score:'),
      hmmpair = keep(i, str_starts, 'hmmpair score:')
    ) %>%
      mutate_all(str_remove, '^.*: ') %>%
      mutate_all(as.numeric) %>%
      mutate(Rfam = j)
  }) %>%
  bind_rows() -> scores.cm.rfam

################################################################################
# collect hmmpair and RNAphylo for the CMfinder detect motifs
# (both biological and random)

cm.motifs %>%
  Sys.glob() %>%
  future_map(function(i) {
    # i <- 'data/H_scores/cmscore/D_search-seqs/K00008_downstream.fna.motif.h1_3.txt'
    # a bit trickier to parse than for the Rfam
    x <- read_lines(i) %>%
      discard(str_detect, ':')
    tibble(
      # extract second word
      sec = x %>%
        strsplit(' ') %>%
        map(2) %>%
        unlist(),
      # extract number part
      num = x %>%
        str_remove('^[^0-9]*') %>%
        # also, for hmmScore remove the second number after ','
        str_remove(',.*$') %>%
        as.numeric
    ) %>%
      spread(sec, num) -> foo
    # ScoreMotif.pl defaults to RNA total, or pair posterior if NA
    if(!('RNA' %in% colnames(foo))) {
      print(i)
      foo$RNA <- foo$pair
    }
    transmute(
      foo,
      RNAphylo = RNA,
      hmmpair = pairSumPf,
      # motif names
      motif = i %>%
        basename %>%
        fs::path_ext_remove(),
      dir = i %>%
        dirname %>%
        basename
    )
  }) %>%
  bind_rows() -> scores.cm.motifs
  

################################################################################
# Load R-scape results

rscape %>%
  Sys.glob() -> rscape.paths 
rscape.paths %>%
  future_map(safely(function(i) {
    x <- read_lines(i)
    
    # Parse length, avg id, etc from MSA line
    x %>%
      keep(str_detect, '^# MSA') %>%
      unlist -> msa
    assertthat::are_equal(length(msa), 1)
    msa %>%
      # remove start with motif name
      str_remove('^# MSA .* (?=nseq)') %>%
      # remove redundant info in paranthesis
      str_remove_all('\\([0-9.]*\\)') %>%
      # conveniently use doulbe space to split key-value pairs
      strsplit('  ') %>%
      unlist %>%
      strsplit(' ') %>%
      map(~ set_names(.x[[2]], .x[[1]])) %>%
      unlist %>%
      as.list() %>%
      as_tibble() -> res
      
    # Attempt to extract observed/expected  (might not exist if zero)
    tryCatch({
      x %>%
        keep(str_detect, '# BPAIRS expected to covary') %>%
        str_remove('^# [ A-Za-z]*') %>%
        str_remove(' \\+/-.*$') %>%
        as.numeric -> res$expected
      
      x %>%
        keep(str_detect, '# BPAIRS observed to covary') %>%
        str_remove('^# [ A-Za-z]*') %>%
        as.numeric -> res$observed
    }, error = function(e) {})
    
    # remember motif and folder
    i %>%
      dirname %>%
      basename -> res$motif
    i %>%
      dirname %>%
      dirname %>%
      basename -> res$dir
    
    return(res)
  })) -> rscape.res

# Check on which parse jobs failed
rscape.res %>%
  map('error') %>%
  map(is.null) %>%
  unlist -> worked
table(worked)

rscape.res %>%
  map('result') %>%
  bind_rows() %>%
  mutate_if(is.numeric, replace_na, 0) -> rscape.res

################################################################################
################################################################################
# Put all together

scores.cm.rfam %>%
  rename(motif = Rfam) %>%
  mutate(dir = 'G_rfam-bacteria-seeds') %>%
  bind_rows(scores.cm.motifs) %>%
  full_join(
    rscape.res %>%
      mutate_at('dir', str_remove, 'H_symlink_'),
    c('motif', 'dir')
  ) -> scores

# Subsitute in-computable scores by Minimum
scores %>%
  mutate_at('RNAphylo', ~ replace_na(.x, min(.x, na.rm = TRUE))) -> scores
# All scores should now be defined
assertthat::are_equal(
  scores %>%  
    filter(!complete.cases(.)) %>%
    nrow,
  0
)


write_tsv(scores, out.scores)

################################################################################
