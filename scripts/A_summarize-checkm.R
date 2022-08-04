# Purpose:
# Create an overview over the genome qualities

library(tidyverse)

# path.tax <- 'data/A_representatives/taxonomy.tsv'
# path.genes <- 'data/A_representatives/genes.tsv.gz'
# path.checkm <- 'data/A_checkm/*/storage/bin_stats_ext.tsv'

path.tax <- unlist(snakemake@input[['tax']])
path.genes <- unlist(snakemake@input[['genes']]) 
# The flag of rules guarantees that these files exist:
path.checkm <- 'data/A_checkm/*/storage/bin_stats_ext.tsv'

# path.fig <- 'data/A_checkm-summary.png'
# path.tbl <- 'data/A_checkm-summary.tsv'

path.fig <- unlist(snakemake@output[['fig']])
path.tbl <- unlist(snakemake@output[['tbl']])

################################################################################
tax <- read_tsv(path.tax)
genes <- read_tsv(path.genes)

################################################################################
# Load a summary of the CheckM results
xs <- Sys.glob(path.checkm)

xs %>%
  map(function(i) {
    # i <- xs[[1]]
    i %>%
      read_lines() %>%
      # issue 1: starts with an invalid identifier
      str_remove('^genome *') %>%
      # issue 2: jsonline does not handle single quotes
      str_replace_all("'", '"') %>%
      jsonlite::parse_json() -> foo
    
    foo$`tax.bio` <- i %>%
      dirname %>%
      dirname %>%
      basename
    
    bar <- c(
      'tax.bio',
      'Genome size',
      'GC', 'GC std',
      'marker lineage',
      '# markers', '# marker sets',
      'Completeness', 'Contamination'
    )
    as_tibble(foo[bar])
  }) %>%
  bind_rows() -> dat
################################################################################
# Bring together with general genome information

dat %>%
  left_join(count(genes, tax.bio, type), 'tax.bio') %>%
  spread(type, n, fill = '0') %>%
  left_join(select(tax, tax.bio, order, species), 'tax.bio') %>%
  separate(tax.bio, c('Taxonomy ID', 'Bioproject'), sep = '\\.') %>%
  mutate_at(c('order', 'species'), str_remove, '^[0-9]+ ') %>%
  mutate_at('order', replace_na, 'unknown') %>%
  select(
    Order = order,
    TxID = `Taxonomy ID`,
    Species = species,
    Bioproject,
    Size = `Genome size`,
    GC,
    GCstd = `GC std`,
    Lineage = `marker lineage`,
    Markers = `# markers`,
    MarkerSets = `# marker sets`,
    Completeness, Contamination,
    Coding = CDS,
    rRNA,
    tRNA,
    everything()
  ) %>%
  arrange(Order, desc(Completeness)) -> dat2

write_tsv(dat2, path.tbl)


################################################################################
# Some pretty plots

dat2 %>%
  mutate(
    Order = fct_lump_n(Order, 4, other_level = 'Other')
  ) -> dat3

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")
  
dat3 %>%
  ggplot(aes(Size / 1e6, as.integer(Coding), color = Order)) +
  # ggsci::scale_color_jco() +
  scale_color_manual(values = cbbPalette) +
  geom_point(size = 4, alpha = .5) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks() +
  xlab('Genome size [Mbp]')  +
  ylab('No. coding genes') +
  theme_bw(18) -> p1

dat3 %>%
  ggplot(aes(GC, as.integer(Coding), color = Order)) +
  # ggsci::scale_color_jco() +
  scale_color_manual(values = cbbPalette) +
  geom_point(size = 4, alpha = .5) +
  xlab('Avg. GC content') +
  ylab('No.coding genes') +
  theme_bw(18) -> p2

dat3 %>%
  ggplot(aes(Completeness, color = Order)) +
  scale_color_manual(values = cbbPalette) +
  stat_ecdf(size = 2, alpha = 0.7) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('Genome completeness') +
  ylab('Empirical cum. density') +
  theme_bw(18) -> p3

dat3 %>%
  ggplot(aes(Size / 1e6, Completeness, color = Contamination)) +
  scale_color_viridis_c() +
  geom_point(size = 5, alpha = .5) +
  scale_x_log10() +
  annotation_logticks(sides = 'b') +
  xlab('Genome size [Mbp]') +
  ylab('Genome completeness') +
  theme_bw(18) -> p4


cowplot::plot_grid(
  p1, p2, p3, p4,
  labels = 'AUTO',
  label_size = 20
) %>%
  ggsave(filename = path.fig, width = 16, height = 11, dpi = 500)
