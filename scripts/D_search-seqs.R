# Purpose: Write out sequences adjacent to the ortholog genes of interest to
# run CMfinder on.

library(tidyverse)
library(plyranges)
library(Biostrings)

library(conflicted)

conflict_prefer("filter", "dplyr")

#############################################################################################
# Snakemake parameters

path.genes <- 'data/A_representatives/genes.tsv.gz'
path.genomes <- 'data/A_representatives/*/genome.fna.gz'
path.trees <- 'data/C_shrink/*/output.tree'

out.intergenic <- unlist(snakemake@output[['intergenic']])
out.seqs <- unlist(snakemake@output[['seqs']])

############################################################################################
# load input data

# Genes
genes <- read_tsv(path.genes)
genes %>%
  select(seqnames = tax.bio.chr, start, end, strand) %>%
         # tax.bio.gene, tax.bio) %>%
  as_granges() -> genes.range
names(genes.range) <- genes$tax.bio.gene

# Genomes
path.genomes %>%
  Sys.glob() %>%
  map(readDNAStringSet) %>%
  purrr::reduce(.f = c) -> genomes

# Tree genes (in this script only the labels of the tips are needed) %
path.trees %>%
  Sys.glob() %>%
  set_names(. %>% dirname %>% basename) %>%
  map(ape::read.tree)  %>%
  map('tip.label') -> tree.genes

#############################################################################################
# Intergenic regions / un-annotated regions 

# length of chromosomes/plastids
tibble(
  seqnames = names(genomes),
  start = 1, end = width(genomes),
  strand = '*'
) %>%
  as_granges() -> genomes.len

# ustranded / without annotation on either strand
unstranded <- setdiff_ranges(genomes.len, reduce_ranges(genes.range))

# stranded / without annotation on at least one strand
bind_ranges(
  mutate(genomes.len, strand = '+'),
  mutate(genomes.len, strand = '-')
) -> genomes.len.both
genomes.len.both %>%
  setdiff_ranges_directed(reduce_ranges_directed(genes.range)) -> stranded

#############################################################################################
# Plot of intergenic region lengths

bind_rows(
  tibble(
    mod = 'unstranded',
    len = width(unstranded)
  ),
  tibble(
    mod = 'stranded',
    len = width(stranded)
  )
) -> dat
  
#2: ECDF
dat %>%
  ggplot(aes(len, color = mod)) +
  ggsci::scale_color_jama(name = NULL) +
  stat_ecdf(size = 1.5) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  scale_x_log10(breaks = c(1, 10, 20, 50, 100, 500, 1e3, 5e3, 1e4)) +
  geom_vline(xintercept = c(20, 500), color = 'red') +
  xlab('Length intergenic regions') +
  ylab('Empirical cum. density') +
  theme_bw(18) +
  theme(legend.position = 'bottom') -> p1

#3: Density
dat %>%
  ggplot(aes(len, color = mod)) +
  ggsci::scale_color_jama(name = NULL) +
  geom_density(size = 1.5) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  scale_x_log10(breaks = c(1, 10, 20, 50, 100, 500, 1e3, 5e3, 1e4)) +
  geom_vline(xintercept = c(20, 500), color = 'red') +
  xlab('Length intergenic regions') +
  ylab('Density') +
  theme_bw(18) +
  theme(legend.position = 'bottom') -> p2

cowplot::plot_grid(p1, p2, nrow = 2,
                   label_size = 18, labels = 'AUTO')

ggsave(out.intergenic, width = 10, height = 14, scale = 1.2)

#############################################################################################
# Only export gaps of lengths within
# generally use max limit region -> for non perfect matches
# afterwards intersect with all gaps and use remainder
# -> Is easier than manually looking up via up down table
LIMITS <- c(20, 500)

#############################################################################################
# Query regions for export

# Genes that have flanking regions with potential to be exported
tree.genes %>%
  unlist %>%
  unique -> tree.genes.uq

# iqtree renamed some genes, by translating special chars to '_'
names(genes.range) <- str_replace_all(
  names(genes.range),
  '[^A-Za-z0-9_.]',
  '_'
)

# anchoring genes for CMfinder
genes.anchors <- genes.range[tree.genes.uq] %>%
  mutate(gene = names(.))

# The potential flanking regions are +/- 500 bp of the anchor gene
c(
  genes.anchors %>%
    flank_upstream(LIMITS[2])  %>%
    mutate(x = 'upstream'),
  genes.anchors %>%
    flank_downstream(LIMITS[2])  %>%
    mutate(x = 'downstream')
) %>%
  mutate(origin.gene = gene) %>%
  select(-gene) %>%
  mutate(potential.i = 1:plyranges::n()) -> potential.flanking

# The flanking sequence have to be
# 1. stranded intergenic
potential.flanking %>%
  join_overlap_intersect_directed(stranded) %>%
  # 2. within the boundary of the plasmid (here, no wrap around origin)
  join_overlap_intersect_directed(genomes.len.both) %>%
  # 3. with the minimal specified size of 20 bp
  filter(width(.) >= LIMITS[1]) %>%
  select(- potential.i) %>%
  unique %>%
  mutate(candidate.i = 1:plyranges::n()) -> candidate.flanking

# ! Challange: there might be multiple gaps within a general range
# Illustration:
#    <GeneX>   <ShortGene>  <AnchorGene>
#        <--------500bp---->
#            ##           ## (two candidate gaps)
#   Solution: Find those closes to the 5'/3' ends of the anchor genes
bind_ranges(
  join_follow_upstream(genes.anchors, candidate.flanking) %>%
    filter(x == 'upstream', gene == origin.gene),
  join_precede_downstream(genes.anchors, candidate.flanking) %>%
    filter(x == 'downstream', gene == origin.gene)
) %>%
  as_tibble() %>%
  select(gene, x, candidate.i) -> to.export

# assert each gene has at most 1 up and 1 downstream region
# (can be fewer due to filtering)
to.export %>%
  count(gene, x) %>%
  filter(n > 1) %>%
  nrow %>%
  assertthat::are_equal(0)

# the sequences to be exported
candidate.flanking %>%
  filter(candidate.i %in% to.export$candidate.i) %>%
  select(candidate.i) -> final.i

to.export.seq <- BSgenome::getSeq(genomes, final.i) 
# assert getSeq returns in same order
assertthat::are_equal(width(to.export.seq), width(final.i))
# for lookup on candidate.i, but make as character
names(to.export.seq) <- as.character(final.i$candidate.i)

#############################################################################################
# Save search sequences for later use
# Extract sequences from genomes

if(!dir.exists(out.seqs)) {
  dir.create(out.seqs)
}

crossing(ko = names(tree.genes), x = c('upstream', 'downstream')) %>%
  mutate(path = sprintf('%s/%s_%s.fna.gz', out.seqs, ko, x)) %>%
  pmap(function(ko, x2, path) {
    # ko <- 'K00003' ; x2 <- 'downstream'
    to.export %>%
      filter(x == x2) %>%
      filter(gene %in% tree.genes[[ko]]) -> foo
    foo %>%
      pull(candidate.i) %>%
      as.character() -> mask
    res <- to.export.seq[mask]
    names(res) <- foo$gene
    # Save but don't create empty files and have at least 5 sequences to work
    if (length(res) >= 5) {
      writeXStringSet(res, path, compress = TRUE)
    }
  })
  
#############################################################################################
