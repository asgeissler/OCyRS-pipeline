# Purpose:
# Inspect alignment of peptide sequences for overal sequence identity 
# consistency, and potential pattern of overall difference for pathways

library(tidyverse)
library(ape)
library(Biostrings)
library(parallel)

library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("first", "dplyr")

# make sure script output is placed in log file
log <- file(unlist(snakemake@log), open="wt")
sink(log)

# the numer of cores to use
cpus <- unlist(snakemake@threads)
# cpus <- 20

###############################################################################

path.og <- unlist(snakemake@input[['og']])
path.kegg <- unlist(snakemake@input[['kegg']])
# path.og <- 'data/B_OGs.tsv'
# path.kegg <- 'data/A_representatives/kegg.tsv.gz'

out.overlap <- unlist(snakemake@output[['overlap']])
out.assoc <- unlist(snakemake@output[['assoc']])
# out.overlap <- 'path-og.jpeg'
# out.assoc <- 'path-og.tsv'

kegg.url <- 'https://rest.kegg.jp/link/pathway/ko'

################################################################################

og <- read_tsv(path.og)
kegg <- read_tsv(path.kegg)

path.ko <- read_tsv(kegg.url, col_names = c('term', 'path')) %>%
  mutate_all(str_remove, '^.*:')

################################################################################
# Check for overlaps between pathways and OGs

# 1. the pathway info
kegg %>%
  filter(db == 'pathway') %>%
  select(path = term, pathway = title, tax.bio.gene) -> path
path %>%
  # 2. overlapping OGs
  right_join(og, 'tax.bio.gene') %>%
  left_join(count(og, term), 'term') %>%
  rename(term.size = n) %>%
  # 3. the propotion of OG contained in the pathway
  count(path, pathway, term, title, term.size) %>%
  rename(overlap = n) %>%
  mutate(prop = overlap / term.size) -> path.og

path.og %>%
  drop_na() %>%
  left_join(
    mutate(path.ko, note = 'yes'),
    c('path', 'term')
  ) %>%
  mutate_at('note', replace_na, 'no') %>%
  ggplot(aes(note, prop)) +
  geom_violin(aes(fill = note)) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = c(
    '#56B4E9',
    '#E69F00'
  )) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('Curated Pathway-KO association') +
  ylab('Proportion of OG overlapping pathway') +
  theme_bw(18) +
  theme(legend.position = 'hide') -> p1

path.og %>%
  drop_na() %>%
  filter(term.size == overlap) %>%
  select(- overlap, - prop) -> path.og2
    
path.og2  %>%
  group_by(path, pathway) %>%
  summarize(
    terms = n(),
    avg.ts = mean(term.size)
  ) %>%
  ungroup() %>%
  mutate(
    # highlight pathways of either top5 avg size or number for terms
    l = ifelse(
      (rank(- avg.ts) <= 5) | (rank(- terms) <= 5),
      pathway, NA_character_)
  ) %>%
  ggplot(aes(terms, avg.ts, label = l)) +
  geom_point(size = 3) +
  ggrepel::geom_label_repel(size = 5) +
  scale_x_log10() +
  xlab('No. OGs in pathway') +
  ylab('Avg. no. genes per OG') +
  theme_bw(18) -> p2
  
cowplot::plot_grid(p1, p2, labels = 'AUTO', label_size = 18)
ggsave(out.overlap, width = 14, height = 7)


# save pathway - OG association for later use
write_tsv(path.og2, out.assoc)

# for reporting

path.og2 %>%
  select(pathway) %>%
  unique %>%
  nrow -> i1

path.og2 %>%
  select(term) %>%
  unique %>%
  nrow -> i2

print(paste(
  i1, 'pathways are associated with',
  i2, 'OGs.'
))
  
################################################################################
################################################################################

# Compute pairwise sequence identities in the alignments

xs <- Sys.glob('data/B_OGs-aln/*.faa.gz')

worker <- function(path) {
  # Load
  x <- readAAMultipleAlignment(path)
  x <- as.AAbin(x)
  # Pairwise distances
  sid <- dist.gene(x, pairwise.deletion = TRUE, method = 'percentage')
  sids <- sid[upper.tri(sid)]
  
  # dist to similarity
  sids <- (1 - sids) * 100
  
  # prep stat
  tibble(
    og = path %>%
      basename() %>%
      str_remove('.faa.gz$'),
    min = min(sids, na.rm = TRUE),
    avg = mean(sids, na.rm = TRUE),
    sd = sd(sids, na.rm = TRUE),
    max = max(sids, na.rm = TRUE)
  )
}

cl <- makeForkCluster(cpus)
res <- parLapply(cl, xs, worker)
stopCluster(cl)

og.sid <- bind_rows(res)

################################################################################


og.sid %>%
  map(~ .x[upper.tri(.x)] ) -> foo

hist(foo$K00012)
foo %>% map(summary)
foo %>% map(mean)
foo %>% map(sd)

################################################################################
################################################################################