# Similarly to G_combine_search.R collect cmsearch results,
# but for G2* the predicted terminator structres
library(tidyverse)
library(furrr)

out.tsv <- unlist(snakemake@output[['tsv']])
out.fig <- unlist(snakemake@output[['fig']])

cpus <- as.integer(unlist(snakemake@threads))
# cpus <- 16
plan(multisession, workers = cpus) 
################################################################################

# names of columns in file
ns <- c('target name', 'accession', 'query name', 'accession2',
  'mdl', 'mdl from', 'mdl to', 'seq from', 'seq to', 'strand',
  'trunc', 'pass', 'gc', 'bias', 'score',
  'E-value', 'inc', 'description of target')

parse <- function(path) {
  ls <- read_lines(path)
  
  # spaces in second line are field breaks
  pos <- str_locate_all(ls[2], ' ')[[1]]
  assertthat::are_equal(pos[, 'start'], pos[, 'end'])
  
  # Extract per non-comment line the fields
  ls %>%
    discard(str_detect, '^#') %>%
    map(function(i) {
      map2(
        c(1, pos[, 'start']),
        c(pos[, 'start'], nchar(ls[2])),
        ~ str_sub(i, .x, .y)
      ) %>%
        map(str_trim) %>%
        unlist()
    }) %>%
    invoke(.f = rbind) %>%
    as_tibble() %>%
    # set column names as from first line comment
    magrittr::set_colnames(ns) %>%
    # remember folder name is tax_bio / genome id
    mutate(
      tax_bio = path %>%
        dirname %>%
        basename
    )
}

# Parse all files
'data/G2_terminators/*.txt' %>%
  Sys.glob() %>%
  future_map(parse) %>%
  invoke(.f = bind_rows) -> tbl

################################################################################

# more concise names
transmute(
  tbl,
  tax_bio = str_remove(`target name`, '\\.[^.]*$'),
  chr = `target name`,
  start = `seq from`, end = `seq to`, strand,
  score, evalue = `E-value`, gc
) -> tbl2

# save 
write_tsv(tbl2, out.tsv)

################################################################################
# A small figure on the no. terminators

'data/A_checkm/checkm_summary.tsv' %>%
  read_tsv() %>%
  mutate(tax_bio = paste(TxID, Bioproject, sep = '.')) -> checkm

tbl2 %>%
  count(tax_bio, name = 'no.term') %>%
  left_join(checkm, 'tax_bio') -> dat

dat %>%
  with(cor.test(no.term, Size)) -> cor.test

dat %>%
  ggplot(aes(Size, no.term, color = GC)) +
  geom_point(size = 3, alpha = .7) +
  scale_color_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  annotate('text', x = 2e6, y = 3e3, color = 'blue', size =5,
           label = sprintf('Pearson-cor.: %+.1f', cor.test$estimate)) +
  geom_smooth(method = 'lm', color = 'blue', se = FALSE) +
  annotation_logticks() +
  xlab('Genome size [bp]') +
  ylab('No. predicted terminators') +
  theme_bw(18)

ggsave(out.fig, width = 7, height = 6)

# > dat$no.term %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 44.0   206.5   461.0   673.8   943.0  3615.0 