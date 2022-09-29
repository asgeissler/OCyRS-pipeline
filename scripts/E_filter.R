# Purpose:
# Check for gappy regions in alignments. These should be removed before
# the SISSIz shuffeling

library(tidyverse)
library(furrr)
library(ape)
library(Biostrings)

in.path <- 'data/D_search-seqs-aln/*.fna.gz'
out.dir <- 'data/E_search-filtered/'
out.gappy <- 'data/E_search-gaps.png'

# cpus <- 16 
cpus <- as.integer(unlist(snakemake@threads))
plan(multisession, workers = cpus) 

THRESHOLD <- .9

################################################################################
# Load alignments

alns <- in.path %>%
  Sys.glob() %>%
  set_names(. %>% basename %>% str_remove('.fna.gz$')) %>%
  future_map(readDNAMultipleAlignment)

################################################################################
# steps in which to count proportion of gaps
xs <- seq(0, 1, .01)

alns %>%
  future_map2(names(.), function(aln, i){
    # i <- 'K00010_downstream'
    # aln <- alns[[i]]
    # the proportion of gaps in the columns
    gap.prop <- consensusMatrix(aln, as.prob = TRUE)['-', ]
    # proportion of columns with prop. of gaps above a cutoff
    col.prop <- cut(gap.prop, xs, include.lowest = TRUE) %>%
      table %>%
      cumsum %>%
      `/`(ncol(aln))
    tibble(
      region = i,
      gappy.threshold = xs[-1] * 100,
      col.prop = col.prop * 100
    )
  }) %>%
  bind_rows() -> gappy

gappy %>%
  ggplot(aes(gappy.threshold, col.prop, group = region)) +
  geom_line(alpha = 0.2, color = 'grey') +
  geom_smooth(aes(group = 1), color = 'blue') +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  geom_vline(xintercept = 90, color = 'red') +
  ylab('Proportion of alignment columnts %') +
  xlab('Threshold % gap positions per column') +
  theme_bw(18)

ggsave(out.gappy, width = 8, height = 6)


################################################################################
# Filter out columns with THRESHOLD or more % gaps

alns %>%
  future_map(function(aln) {
    aln %>%
      maskGaps(min.fraction = THRESHOLD) %>%
      as("DNAStringSet") -> aln2
    # remeber to also remove sequences that might now only consist of gaps
    mask <- letterFrequency(aln2, '-') / nchar(aln2) < .999
    aln2[names(aln2)[mask]]
  }) -> alns.filtered

# Exclude now empty regions
alns.filtered <- alns.filtered %>%
  keep(~ length(.x) > 0)

# Save
if(!dir.exists(out.dir)) {
  dir.create(out.dir)
}
alns.filtered %>%
  future_map2(names(.), ~ writeXStringSet(
    .x,
    sprintf('%s/%s.fna.gz', out.dir, .y),
    compress = TRUE
  ))

################################################################################
