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

cpus <- 10
plan(multicore, workers = cpus)

################################################################################

in.rfam <- 'data/G_rfam-cmsearch.tsv.gz'
in.term <- 'data/G2_terminators.tsv.gz'
in.motifs <- 'data/J_motif-aln-seq-pos.tsv'
in.cats = 'data/J_FDR-categories.tsv'

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

write_tsv(search.tbl, 'foo.tsv.gz')

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
  mutate(region = search.ranges2$region)

ranges2 %>%
  as_tibble() %>%
  select(- ranges.row) %>%
  write_tsv('foo.tsv.gz')

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

ggsave('foo2.jpeg',
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
  )
  mutate(
    recall = no.recall / no.pos,
    precision = no.recall / no.motif,
    # F1 = 2 TP / (2* TP + FP + FN)
    F1 = 2 * no.recall / ( 2 * no.recall + (no.seqs.search - no.recall) + (motif.seqs - no.recall))
  ) -> over.stat
over.stat %>%
  ggplot(aes(recall, precision, color = type2)) +
  # ggsci::scale_color_jco(name = NULL) +
  scale_color_manual(values = cbbPalette, name = NULL) +
  geom_point(size = 3, alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  xlab('recall (overlapping motifs rel. to\nno. RNIE/Rfam model hits)') +
  ylab('precision (overlapping motifs\nrel. to no. motif homologs') +
  annotation_logticks() +
  theme_bw(18) +
  theme(legend.position = 'bottom') -> p

p <- ggExtra::ggMarginal(p, type = 'histogram')
ggsave('foo3.jpeg', plot = p, width = 10, height = 10, dpi = 500)

over.stat %>%
  group_by(motif) %>%
  slice_max(recall) %>%
  left_join(cats, 'motif') %>%
  ggplot(aes(recall, precision, color = category)) +
  geom_point() + facet_wrap(~ type2, scales = 'free')
  ggplot(aes(type2, precision, color = category)) +
  geom_boxplot()
categroy

###

ns <- sprintf('Overlapping RNA structure motifs (≥ %g bp)', cutoff.recall.overlap)
over.stat %>%
  select(motif, type2) %>%
  unique %>%
  count(type2) %>%
  arrange(desc(n)) %>%
  rename(!! ns:= n) %>%
  knitr::kable()

################################################################################
# Select the potential novel homologs


no.seq %>%
  filter(type == 'CMfinder motif') %>%
  select(- type) %>%
  rename(motif = name) %>%
  anti_join(over.stat, 'motif') -> potential.novel


ranges.nonredundant %>%
  filter(type == 'CMfinder motif') %>%
  filter(name %in% potential.novel$motif) %>%
  select(- type, - rfam) -> potential.novel.ranges

################################################################################
# Compare KEGG terms, associated pathways when comparing only the
# alignment sequences or also additional homologs

# Questions of interest:
# 1. Pathways associated via KO number
# 2. Pathways associated for genes +/- 500 bp away to include sequence homologs
# 3. When encluding cmsearch predicted homologs

# Q1
'https://rest.kegg.jp/link/pathway/ko' %>%
  read_tsv(col_names = c('ko', 'path')) %>%
  filter(str_detect(path, 'map')) %>%
  mutate_at('ko', str_remove, '^ko:') %>%
  mutate_at('path', str_remove, '^path:') -> ko.path

'https://rest.kegg.jp/list/pathway' %>%
  read_tsv(col_names = c('path', 'path_name')) %>%
  filter(str_detect(path, 'map')) %>%
  mutate_at('path', str_remove, '^path:') -> path.names

potential.novel %>%
  mutate(ko = str_remove(motif, '_.*$')) %>%
  left_join(ko.path, 'ko') %>%
  left_join(path.names, 'path') %>%
  select(path, path_name, motif) %>%
  unique %>%
  drop_na %>%
  count(
    path, path_name,
    name = 'Motifs associated via CMfinder search anchor'
  ) -> q1


# Q2

'data/A_representatives/kegg.tsv.gz' %>%
  read_tsv() -> gene.kegg

'data/A_representatives/genes.tsv.gz' %>%
  read_tsv() -> genes

# search in +/- 500
genes %>%
  select(
    seqnames = tax.bio.chr,
    start, end, strand,
    tax.bio.gene
  ) %>%
  plyranges::as_granges() %>%
  mutate(
    start = start - 500,
    end = end + 500
  ) -> genes500
  
# check overlaps
genes500 %>%
  plyranges::join_overlap_intersect(potential.novel.ranges) %>%
  as_tibble() -> genes500.novel.over


genes500.novel.over %>%
  left_join(gene.kegg, 'tax.bio.gene') %>%
  filter(db == 'pathway') %>%
  select(
    motif = name, tax.bio.gene, alignment.seq,
    path = term, path_name = title
  ) %>%
  unique -> genes500.motif.path
genes500.motif.path %>%
  filter(alignment.seq) %>%
  drop_na %>%
  count(
    path, path_name,
    name = 'Motifs with align. seq. near genes'
  ) -> q2

# Q3, repeated but without the alignment.seq filter
genes500.motif.path %>%
  drop_na %>%
  count(
    path, path_name,
    name = 'Motifs when incl. cmsearch homologs'
  ) -> q3


list(q1, q2, q3) %>%
  reduce(full_join, by = c('path', 'path_name')) %>%
  View

### Why does the number decrease for map00543?
potential.novel %>%
  mutate(ko = str_remove(motif, '_.*$')) %>%
  left_join(ko.path, 'ko') %>%
  left_join(path.names, 'path') %>%
  filter(path == 'map00543')

gene.kegg %>%
  filter(term == 'K19294') %>%
  select(tax.bio.gene) %>%
  left_join(gene.kegg) %>%
  filter(term == 'map00543')
# Problem: KEGG rest api lists KO-path associations not in proGenomes :*(

################################################################################

potential.novel %>%
  left_join(homologs, 'motif') %>%
  filter(str_detect(tax.bio.chr, '^32049')) %>%
  select(motif, tax.bio.chr, start, end, strand, alignment.seq)

potential.novel %>%
  left_join(homologs, 'motif') %>%
  filter(str_detect(tax.bio.chr, '^32049')) %>%
  select(motif) %>%
  left_join(fdr, 'motif')


################################################################################

