# Purpose:
# Check for overlaps between motif alignment sequence positions and cmsearch
# hits of
# (A) bacterial Rfam
# (B) RNIE terminator predictions

library(tidyverse)
library(plyranges)
library(furrr)

library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("n", "dplyr")
conflict_prefer("first", "dplyr")


# cpus <- as.integer(unlist(snakemake@threads))
cpus <- 10
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
# (Omitted to explicitly save and it is now easier to re-implement than risking
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
  future_map(Biostrings::readDNAStringSet) %>%
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
# Inspect Evalues to estimate error of CMsearch ahead of CMfinder recall

list(ranges, ranges2) %>%
  map(function(i) {
    i %>%
      filter(!is.na(evalue)) %>%
      as_tibble() %>%
      mutate(
        name = replace_na(name, 'RNIE'),
        rfam = replace_na(rfam, 'Terminator')
      ) %>%
      group_by(name, rfam) %>%
      summarize(
        no.hits = n(),
        evalue.sum = sum(evalue),
        cmsearch.fdr = evalue.sum / no.hits
      )
  }) %>%
  invoke(.f = left_join, by = c('name', 'rfam')) %>%
  arrange(cmsearch.fdr.y) %>%
  transmute(
    name, rfam,
    #
    'CMsearch hits' = as.integer(no.hits.x),
    'Sum of E-Values' = evalue.sum.x,
    'FDR' = cmsearch.fdr.x,
    #
    'Hits in search regions' = as.integer(replace_na(no.hits.y, 0)),
    'Corresponding E-Values' = evalue.sum.y,
    'Expected FDR for CMfinder' = cmsearch.fdr.y
  ) -> rfam.fdr

################################################################################
# Compute overlaps, but only per group to prevent erroneous multiple counts

ranges2 %>%
  filter(type == 'Candidate motif') %>%
  reduce_ranges_directed() %>%
  mutate(reduced.row = 1:plyranges::n()) -> ranges2.red.motifs

ranges2 %>%
  filter(type != 'Candidate motif') %>%
  group_by(type, rfam) %>%
  reduce_ranges_directed() %>%
  ungroup %>%
  mutate(reduced.row = 1:plyranges::n()) -> ranges2.red.ref

join_overlap_intersect(
  ranges2.red.motifs,
  ranges2.red.ref %>%
    mutate(strand2 = strand)
) %>%
  mutate(
    overlap = IRanges::width(.),
    orientation = ifelse(
      strand == strand2,
      'sense',
      'anti-sense'
    )
  ) -> overlaps

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
# Requested: "Recall" if only a single position overlaps...

overlaps %>%
  as_tibble() %>%
  select(rfam) %>%
  unique %>%
  mutate(
    rfam = replace_na(rfam, 'Terminator'),
    'Overlaps motif' = 'yes'
  ) %>%
  right_join(rfam.fdr, 'rfam') %>%
  left_join(
    rfam.types %>%
      select(name = Family, rfam.type = type),
    'name'
  ) -> recall.single.overlap

################################################################################
# Improved: Recall relative to actual positions

list(
  ranges2.red.ref %>%
    as_tibble() %>%
    mutate(rfam = replace_na(rfam, 'Terminator')) %>%
    select(rfam, reduced.row) %>%
    unique %>%
    count(rfam, name = 'reduced.ranges.ref'),
  overlaps %>%
    as_tibble() %>%
    filter(orientation == 'sense') %>%
    mutate(rfam = replace_na(rfam, 'Terminator')) %>%
    select(rfam, reduced.row.y) %>%
    unique %>%
    count(rfam, name = 'with.sense.overlap'),
  overlaps %>%
    as_tibble() %>%
    filter(orientation == 'anti-sense') %>%
    mutate(rfam = replace_na(rfam, 'Terminator')) %>%
    select(rfam, reduced.row.y) %>%
    unique %>%
    count(rfam, name = 'with.anti.overlap')
) %>%
  purrr::reduce(.f = full_join, by = 'rfam') %>%
  mutate_if(is.numeric, replace_na, 0) %>%
  mutate(
    sense.recall = with.sense.overlap / reduced.ranges.ref,
    anti.recall = with.anti.overlap / reduced.ranges.ref
  ) -> over.stat

recall.single.overlap %>%
  left_join(over.stat, 'rfam') %>% 
  transmute(
    Family = name,
    Name = rfam,
    Type = replace_na(rfam.type, 'Terminator'),
    `CMsearch hits`, `Sum of E-Values`, `FDR`,
    `Hits in search regions`, `Corresponding E-Values`,
    `Expected FDR for CMfinder`,
    `Overlaps motif`,
    'CMfinder sense recall' = sense.recall,
    'Anti-sense recall' = anti.recall
  ) -> recall.full

recall.full %>%
  write_tsv('data/J_novel/reference-motif-overlap-stats.tsv')

################################################################################
# Plot FDR expected for CMsearch

recall.full %>%
  select(Name, Type, 'CMsearch FDR, overall' = FDR,
         'FDR within search regions' = `Expected FDR for CMfinder`) %>%
  drop_na %>%
  gather('k', 'v', - c(Name, Type)) %>%
  mutate(Type = fct_reorder(Type, v)) %>%
  ggplot(aes(Type, v, color = k, group = paste(Type, k))) +
  geom_boxplot(position = position_dodge(preserve = 'single')) +
  geom_hline(color = 'red', yintercept = 1e-3) +
  annotate('text',
           x = 0.5, y = 0.2, hjust = 0,
           label = 'FDR = 0.001',
           color = 'red', size = 5) +
  ggsci::scale_color_jama(name = NULL) +
  xlab('RNA family') +
  ylab('FDR estimate\nSum of E-values over no. CMsearch hits') +
  scale_y_log10() +
  guides(color = guide_legend(direction = 'vertical')) +
  theme_bw(18) +
  theme(legend.position = 'bottom') -> p1



################################################################################
# Recall stats

recall.full %>%
  filter(`Overlaps motif` == 'yes') %>%
  gather('recall', 'v', c(`CMfinder sense recall`, `Anti-sense recall`)) %>%
  mutate(Type = fct_reorder(Type, v, .fun = max)) %>%
  ggplot(aes(Type, v, color = fct_rev(recall))) +
  geom_boxplot() +
  xlab('RNA family') +
  ylab('Recall of CMfinder strucutres positions\nrelative to CMsearch hits in search regions') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  scale_color_manual(values = c('blue', 'red'), name = NULL) +
  # ggsci::scale_color_d3(name = NULL) +
  guides(color = guide_legend(direction = 'vertical')) +
  theme_bw(18) +
  theme(legend.position = 'bottom') -> p2

cowplot::plot_grid(
  p1, p2, 
  labels = 'AUTO', label_size = 16
)

ggsave('data/J_novel/recall-barplots.jpg', width = 17, height = 8, dpi = 400)

################################################################################
# Overview table on recall

list(
  recall.full %>%
    count(Type, name = 'RNA families with CMsearch hits in Cyanobacteria'),
  recall.full %>%
    filter(`Hits in search regions` > 0) %>%
    count(Type, name = 'With hits in search regions'),
  recall.full %>%
    filter(`Overlaps motif` == 'yes') %>%
    count(Type, name = 'Overlapping a CMfinder motif'),
  recall.full %>%
    filter(`CMfinder sense recall` >= 0.5) %>%
    count(Type, name = 'CMfinder recalls ≥50% of CMsearch hits, same strand'),
  recall.full %>%
    filter(`Anti-sense recall` >= 0.5) %>%
    count(Type, name = 'Anti-sense recall ≥50%')
) -> foo

foo %>%
  purrr::reduce(.f = full_join, 'Type') %>%
  gather('k', 'v', - Type) %>%
  spread(Type, v) %>%
  mutate_at('k', fct_relevel, foo %>%
              map(colnames) %>%
              unlist %>%
              discard(~ .x == 'Type') %>%
              unlist) %>%
  mutate_if(is.numeric, replace_na, 0) %>%
  arrange(k) %>%
  select(
    'Description' = k,
    sRNA, `Cis-reg`,
    other, CRISPR,
    rRNA, tRNA,
    everything()
  ) %>%
  write_tsv('data/J_novel/overview.tsv')

################################################################################
# Motif-centric overview table of
# - Predominant recalling family type per motif category
# - Novel motifs
# - Motifs in search regions without references (or at least Rfam)

# Motifs overlapping a reduce motif that overlaps a reference
ranges2 %>%
  filter(type == 'Candidate motif') %>%
  join_overlap_intersect(
    # Reduced motif ranges that overlap a reference
    ranges2.red.motifs[unique(overlaps$reduced.row.x)]
  )  -> mot.red.ref
# Motifs that don't overlap references, or reduced motifs that overlapped
ranges2 %>%
  filter(type == 'Candidate motif') %>%
  filter(! (name %in% mot.red.ref$name) ) -> pot.novel
# Potential novel hits per search regions with / without overlapping 
# CMsearch hits
{
  # List of intergenic regions with references
  search.ranges2 %>%
    filter(type != 'Candidate motif') %>%
    as_tibble() %>%
    pull(intergenic.region) %>%
    unique -> inter.ref
  # motifs with ≥1 hit in such regions
  search.ranges2 %>%
    filter(type == 'Candidate motif') %>%
    filter(intergenic.region %in% inter.ref) %>%
    as_tibble() %>%
    pull(name) %>%
    unique -> inter.ref.mot
  # Negative search, potential novel motifs without such hits
  pot.novel %>%
    mutate(search.region = ifelse(
      name %in% inter.ref.mot, 
      'has CMsearch ref',
      'w/o ref'
    ))
} -> pot.novel.status

list(
  cats %>%
    count(category, name = 'No. CMfinder motifs'),
  cats %>%
    filter(motif %in% mot.red.ref$name) %>%
    count(category, name = 'Motifs recalling a reference'),
  cats %>%
    filter(motif %in% pot.novel$name) %>%
    count(category, name = 'Motifs without reference overlaps'),
  cats %>%
    filter(motif %in%  unlist(
      pot.novel.status[pot.novel.status$search.region == 'w/o ref']$name
    )) %>%
    count(category, name = 'In search regions without CMsearchhits')
) %>%
  purrr::reduce(.f = full_join, by = 'category') -> foo


foo %>%
  select_if(is.numeric) %>%
  map(sum) %>%
  c('category' = 'Total') %>%
  as_tibble %>%
  bind_rows(foo, .) %>%
  write_tsv('data/J_novel/overview-motifs.tsv')

################################################################################
# Select the potential novel motifs

pot.novel.status %>%
  as_tibble() %>%
  select(motif = name, search.region) %>%
  unique %>%
  left_join(cats, 'motif') -> potential.novel

write_tsv(potential.novel, 'data/J_novel/potentially-novel-motifs.tsv')

################################################################################

