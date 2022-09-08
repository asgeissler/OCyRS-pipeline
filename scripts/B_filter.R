# Prefilter aligned sequences
# Remove those sequences that have predominately all amino acids
# in regions that are gappy
library(tidyverse)
library(Biostrings)

path <- unlist(snakemake@input)

path.out <- unlist(snakemake@output[['res']])
path.fig <- unlist(snakemake@output[['fig']])

# make sure script output is placed in log file
log <- file(unlist(snakemake@log), open="wt")
sink(log)

################################################################################

seqs <- readAAStringSet(path)
aln <- AAMultipleAlignment(seqs)

################################################################################
# mask clumns with 80% gaps
aln.masked <- maskGaps(aln, min.fraction = .8)

################################################################################
# Compare proportion of AAs in/outside of masked blocks

aas <- alphabetFrequency(aln)[, AA_STANDARD] %>%
  apply(1, sum)
aas.in.blocks <- alphabetFrequency(aln.masked)[, AA_STANDARD] %>%
  apply(1, sum)

out.ratio <- 1 - aas.in.blocks / aas

################################################################################
# remove outlier
mask <- out.ratio > .2

# clean columns that now only have gaps
seqs[!mask] %>%
  AAMultipleAlignment() %>%
  maskGaps(min.fraction = 1) %>%
  as("AAStringSet") -> res

res.aln <- AAMultipleAlignment(res)

################################################################################
### Compile plot for future debugging
tibble(x = out.ratio) %>%
  ggplot(aes(x)) +
  # geom_histogram() +
  stat_ecdf() +
  xlab('Proportion of AAs in masked columns') +
  geom_vline(xintercept = .1, color = 'red') +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  theme_bw(18) +
  ggtitle(
    path,
    sprintf(
      paste(
        'Alignment had %s sequences and %s columns',
        '%s columns consist of >80%% gaps',
        '%s sequences have >90%% of AAs in gapped regions',
        sep = '\n'
      ),
      nrow(aln), ncol(aln),
      maskedncol(aln.masked),
      sum(mask)
      
    )
  ) -> p

################################################################################
# Save results
writeXStringSet(res, path.out, compress = TRUE)
ggsave(path.fig, p, width = 9, height = 8)

