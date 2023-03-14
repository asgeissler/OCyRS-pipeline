# Purpose:
# Check for overlaps between motif alignment sequence positions and cmsearch
# hits of
# (A) bacterial Rfam
# (B) RNIE terminator predictions

library(tidyverse)
library(tidygraph)
library(plyranges)
library(furrr)

library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("n", "dplyr")
conflict_prefer("first", "dplyr")


cpus <- as.integer(unlist(snakemake@threads))
# cpus <- 10
plan(multicore, workers = cpus)

################################################################################

in.rfam <- 'data/G_rfam-cmsearch.tsv.gz'
in.term <- 'data/G2_terminators.tsv.gz'
in.motifs <- 'data/J_motif-aln-seq-pos.tsv'
in.cats <- 'data/J_FDR-categories.tsv'

path.genes <- 'data/A_representatives/genes.tsv.gz'
path.genomes <- 'data/A_representatives/*/genome.fna.gz'
path.trees <- 'data/C_shrink/*/output.tree'

# Colorblind-friendly palettes of the Color Universal Design
# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")

################################################################################
# Re-create search sequence coordinates computation
# (Omitted to explicilty save and it is now easier to re-implement than risking
#  a broken Snakemake pipeline)
# Code copy/pasted from `D_search-seqs.R` (see comment in that file for explanations)

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
  map(Biostrings::readDNAStringSet) %>%
  purrr::reduce(.f = c) -> genomes

# Tree genes (in this script only the labels of the tips are needed) %
path.trees %>%
  Sys.glob() %>%
  set_names(. %>% dirname %>% basename) %>%
  future_map(function(i) {
    ape::read.tree(i)$tip.label
  }) -> tree.genes
  # map(ape::read.tree)  %>%
  # map('tip.label') -> tree.genes

# match genes to KO terms
tree.genes %>%
  map2(names(.), ~ tibble(gene = .x, KO = .y)) %>%
  bind_rows() -> ko.gene

# length of chromosomes/plastids
tibble(
  seqnames = names(genomes),
  start = 1, end = width(genomes),
  strand = '*'
) %>%
  as_granges() -> genomes.len

# stranded / without annotation on at least one strand
bind_ranges(
  mutate(genomes.len, strand = '+'),
  mutate(genomes.len, strand = '-')
) -> genomes.len.both
genomes.len.both %>%
  setdiff_ranges_directed(reduce_ranges_directed(genes.range)) -> stranded

LIMITS <- c(20, 500)

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
  select(gene, x, candidate.i) %>%
  # Divergend from D_*R script from here onward
  left_join(ko.gene, 'gene') %>%
  unite(region, KO, x) -> to.intergenic

candidate.flanking[to.intergenic$candidate.i] %>%
  mutate(region = to.intergenic$region) %>%
  select(region, origin.gene, intergenic.region = candidate.i) -> search.ranges

search.ranges %>%
  as_tibble() -> search.tbl

write_tsv(search.tbl, 'data/J_novel/all_intergenic_regions.tsv.gz')

################################################################################
# Load remaining data

rfam <- read_tsv(in.rfam)
term <- read_tsv(in.term)
motifs <- read_tsv(in.motifs)
cats <- read_tsv(in.cats)


################################################################################
# as genomic range

# Process motifs and RNIE/terminators equally
list(
  'RNIE terminator' = term %>%
    select(- tax_bio, - gc),
  Rfam = rfam %>%
    select(- tax_bio, - gc) %>%
    rename(name = rf, rfam = name),
  'Candidate motif' = motifs %>%
    rename(name = motif)
) %>%
  map2(names(.), ~ mutate(.x, type = .y)) %>%
  bind_rows() %>%
  mutate(ranges.row = 1:n()) %>%
  rename(seqnames = chr) %>%
  mutate(
    tmp = start,
    start = ifelse(start < end, start, end),
    end = ifelse(tmp < end, end, tmp)
  ) %>%
  select(- tmp) %>%
  plyranges::as_granges() %>%
  mutate(len = IRanges::width(.)) -> ranges


################################################################################
# Limit to search region
# also filter out alignment sequence positions outside of the initial search
# region (downside of the exact sequence lookup)


search.ranges %>%
  join_overlap_inner(ranges) %>%
  # Extract KO info of motifs
  mutate( motif.ko = str_replace(
    name,
    '^(K[0-9]{5}_[updown]{2,4}stream).fna.motif.*$',
    '\\1'
  )) %>%
  # Keep all Rfam / RNIE or motifs (but then only if search region match)
  filter(
    (type != 'Candidate motif') | (motif.ko == region)
  ) -> search.ranges2

ranges2 <- ranges[search.ranges2$ranges.row] %>%
  mutate(region = search.ranges2$region) %>%
  unique()

ranges2 %>%
  as_tibble() %>%
  select(- ranges.row) %>%
  write_tsv('data/J_novel/references_inside-of_intergenic_regions.tsv.gz')

################################################################################
# Compute overlaps, but only per group to prevent erroneous multiple counts

join_overlap_intersect(
  ranges2 %>%
    filter(type == 'Candidate motif') %>%
    select(region, motif = name, ranges.row.motif = ranges.row),
  ranges2 %>%
    filter(type != 'Candidate motif') %>%
    mutate(strand2 = strand) %>%
    select(- len)
) %>%
  mutate(
    overlap = IRanges::width(.),
    orientation = ifelse(
      strand == strand2,
      'sense',
      'anti-sense'
    )
  ) %>%
  filter(region.x == region.y) %>%
  select(
    region = region.x, motif, ranges.row.motif,
    overlap, orientation,
    type, name, rfam, ranges.row, score, evalue
  ) -> overlaps

################################################################################
# Determine bp cutoff for recall

cutoff.recall.overlap <- 5

overlaps %>%
  as_tibble %>%
  ggplot(aes(overlap, color = orientation)) +
  stat_ecdf(size = 1.2) +
  ggsci::scale_color_jama() +
  # scale_x_continuous(breaks = c(seq(0, 50, 10), 100, 150)) +
  scale_x_log10(breaks = c(1, 5, 10, 50, 100, 150)) +
  annotation_logticks(sides = 'b') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_vline(xintercept = cutoff.recall.overlap, color = 'blue') +
  xlab('Overlap with candidate motif [bp]') +
  ylab('Empirical cumulative density') +
  facet_wrap( ~ type) +
  theme_bw(18) +
  theme(legend.position = 'bottom')

ggsave('data/J_novel/reference-motif-overlaps.jpg',
       width = 12, height = 6, dpi = 500)

################################################################################
# Ahead of recall, query Rfam family types for more informative plots

'https://raw.githubusercontent.com/Rfam/rfam-taxonomy/master/domains/bacteria.csv' %>%
  read_csv() %>%
  select(Family, Name = `Rfam ID`, type = `RNA type`) -> rfam.types

keep.type <- c('Cis-reg', 'antisense', 'CRISPR', 'ribozyme', 'rRNA', 'tRNA', 'sRNA')

rfam.types %>%
  mutate(
    type2 = type %>%
      map(function(i) {
        keep.type %>%
          map(~ str_detect(i, .x)) %>%
          unlist %>%
          which %>%
          first -> i
        keep.type[i]
      }) %>%
      unlist
  ) %>%
  mutate_at('type2', replace_na, 'other') %>%
  rename(type.full = type, type = type2) -> rfam.types

################################################################################
# Determine recall, but specific per search regions!

ranges2 %>%
  as_tibble() %>%
  count(region, type, name, name = 'no.pos') -> no.pos


overlaps %>%
  filter(overlap >= cutoff.recall.overlap) %>%
  as_tibble %>%
  left_join(rename(rfam.types, rf.type = type), c('name' = 'Family')) %>%
  mutate(
    type2 = case_when(
      type == 'RNIE terminator' ~ 'RNIE terminator',
      TRUE ~ paste('Rfam', rf.type)
    )
  ) %>%
  count(
    region, motif,
    orientation, type, type2, name,
    name = 'no.recall'
  ) %>%
  left_join(no.pos, c('region', 'type', 'name')) %>%
  left_join(
    select(no.pos, motif = name, region, no.motif = no.pos),
    c('region', 'motif')
  ) %>%
  mutate(
    recall = no.recall / no.pos,
    precision = no.recall / no.motif,
  ) -> over.stat

over.stat %>%
  write_tsv('data/J_novel/reference-motif-overlap-stats.tsv')

################################################################################

over.stat %>%
  select(motif, orientation, type2, name, recall) %>%
  left_join(cats, 'motif') %>%
  mutate_at('category', str_replace, ' \\(', '\n(') %>%
  mutate(type3 = fct_reorder(type2, recall, .desc = TRUE)) %>%
  ggplot(aes(type3, recall, color = orientation)) +
  # scale_color_manual(values = cbbPalette, name = NULL) +
  ggsci::scale_color_jama() +
  geom_boxplot(position = position_dodge2(preserve = 'single')) +
  xlab(NULL) +
  ylab('Recall rates per predicted RNA structure motifs\nconstrained to ortholog search regions') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  facet_grid(~ category, scales = 'free_x', space = 'free_x') +
  theme_bw(18) +
  theme(
    legend.position = 'bottom',
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) -> p1


################################################################################

over.stat %>%
  left_join(cats, 'motif') %>%
  mutate_at('category', str_replace, ' \\(', '\n(') %>%
  mutate(precision = pmin(precision, 1)) %>%
  ggplot(aes(recall, precision, color = type2)) +
  geom_point(size = 3, alpha = 0.7) +
  ggsci::scale_color_igv(name = NULL) +
  facet_grid(orientation ~ category) +
  theme_bw(18) +
  theme(legend.position = 'bottom') -> p2

cowplot::plot_grid(p1, p2, labels = 'AUTO', label_size = 16, cols = 1)

ggsave('data/J_novel/recall-precision-plot.jpg', width = 16, height = 18, dpi = 400)

################################################################################

over.stat %>%
  left_join(cats, 'motif') %>%
  filter(recall > .5, precision > .5) %>%
  # select(type2, name, orientation, motif, category) %>%
  select(type2, name, motif, category) %>%
  unique -> foo

full_join(
  foo %>%
    select(- motif) %>%
    unique %>%
    count(type2, category, name = 'Families'),
  foo %>%
    select(- name) %>%
    unique %>%
    count(type2, category, name = 'Motifs'),
  c('type2', 'category')
) %>%
  unite('nice', Families, Motifs, sep = ' - ') %>%
  spread(category, nice, fill = '0 - 0') %>%
  rename('Recall & precision > 0.5; No. families - motifs' = type2) %>%
  write_tsv('data/J_novel/overview.tsv')
  

################################################################################
# Select the potential novel motifs

cats %>%
  anti_join(over.stat, 'motif') -> potential.novel

write_tsv(potential.novel, 'data/J_novel/potentially-novel-motifs.tsv')

################################################################################

