# Purpose:
# Look up the coordinates of the detect motifs in the bacterial genomes

library(tidyverse)
library(Biostrings)
library(plyranges)
library(furrr)

in.path <- 'data/I_fdr.tsv'
in.genes <- 'data/A_representatives/genes.tsv.gz'
in.rfam <- 'data/G_rfam-cmsearch.tsv.gz'

# make sure script output is placed in log file
log <- file(unlist(snakemake@log), open="wt")
sink(log)

# the numer of cores to use
# cpus <- 16
cpus <- as.integer(unlist(snakemake@threads))
plan(multisession, workers = cpus) 

################################################################################

# paths to motifs
motifs <- in.path %>%
  read_tsv() %>%
  select(motif, region) %>%
  separate(region, c('og', 'side'), sep = '_', remove = FALSE) %>%
  mutate(path = sprintf('data/F_cmfinder/D_search-seqs/%s/%s', region, motif))

# gene coordinetes
genes <- read_tsv(in.genes) %>%
  select(-type, -gc.content, - product, - gene) %>%
  dplyr::rename(gene.start = start, gene.end = end, gene.strand = strand) %>%
  # iqtree renamed some genes, by translating special chars to '_'
  mutate_at('tax.bio.gene', str_replace_all, '[^A-Za-z0-9_.]', '_')

################################################################################

# Nuisance: Determining positions for upstream cases (see pos.helper comments)
#           requires knowledge of the full length of search sequences.
# -> Load lengths

reglen.help <- function(region, path) {
  xs <- readDNAStringSet(path)
  tibble(
    tax.bio.gene = names(xs),
    region.len = width(xs)
  ) %>%
    mutate_at('tax.bio.gene', str_replace_all, '[^A-Za-z0-9_.]', '_') %>%
    mutate(region = region)
}

region.lengths <- motifs %>%
  select(region) %>%
  unique %>%
  mutate(path = sprintf('data/D_search-seqs/%s.fna.gz', region)) %>%
  future_pmap(reglen.help) %>%
  bind_rows()

################################################################################
# Illustration on Up/downstream location relative to +/- located genes.
#
# Fwd:    Upstream      GeneA  Downstream
# Genome: ======================================================================
# Rev:                                        Downstream    GeneB   Upstream

# Assumption: CMfinder's relative positions are one-based and inclusive
#
# Helper implementation to compute absolute positions

pos.helper <- function(x) {
  # x is tibble with keys:
  # - gene.start
  # - gene.end
  # - gene.strand (+/- as chr)
  # - rel.start
  # - rel.end
  # - region.len
  # - side is either upstream or downstream as chr
  mutate(
    x, 
    start = case_when(
      (side == 'upstream')   & (gene.strand == '+') ~ gene.start - region.len + rel.start - 1,
      (side == 'upstream')   & (gene.strand == '-') ~ gene.end   + region.len - rel.end + 1,
      (side == 'downstream') & (gene.strand == '+') ~ gene.end                + rel.start,
      (side == 'downstream') & (gene.strand == '-') ~ gene.start              - rel.end
    ),
    end = case_when(
      (side == 'upstream')   & (gene.strand == '+') ~ gene.start - region.len + rel.end - 1,
      (side == 'upstream')   & (gene.strand == '-') ~ gene.end   + region.len - rel.start + 1,
      (side == 'downstream') & (gene.strand == '+') ~ gene.end                + rel.end,
      (side == 'downstream') & (gene.strand == '-') ~ gene.start              - rel.start
    )
  )
  # Check
  # +down OK
  # +up
  # -down
  # -up
}

################################################################################

worker <- function(path, side, region, motif) {
  # path <- 'data/F_cmfinder/D_search-seqs/K00003_upstream/K00003_upstream.fna.motif.h1_1'
  # side <- 'upstream'
  # region <- 'K00003_upstream'
  # Load the relative positions from the Stockholm motif file
  rel <- path %>%
    read_lines %>%
    keep(str_detect, '^#=GS') %>%
    keep(str_detect, ' DE ') %>%
    str_remove('#=GS ') %>%
    str_remove('\t.*$') %>%
    tibble(info = .) %>%
    separate(info, c('tax.bio.gene', 'xy'), sep = ' +DE +') %>%
    separate(xy, c('rel.start', 'rel.end'), sep = '\\.\\.') %>%
    mutate(across(starts_with('rel'), as.numeric))
  
  
  # use helper to compute absolute values
  rel %>%
    left_join(genes, 'tax.bio.gene') %>%
    mutate(region = region, side = side) %>%
    left_join(region.lengths, c('region', 'tax.bio.gene')) %>%
    pos.helper() %>%
    dplyr::rename(strand = gene.strand) %>%
    mutate(motif = motif)
}


motifs %>%
  select(path, side, region, motif) %>%
  future_pmap(worker) %>%
  bind_rows() -> abs.pos


assertthat::are_equal(
  abs.pos %>%
    filter(start > end) %>%
    nrow,
  0
)

################################################################################
################################################################################
# Though not strictly needed, select some upstream/downstream cases
# and double check that the computed position are matching the sequences.
# Doing such a test ensures that the assumption on CMfinder's output is given.

abs.pos %>%
  mutate(side = ifelse(str_detect(motif, '_upstream.fna.motif'),
                       'upstream', 'downstream')) %>%
  # group_by(side, strand) %>%
  # slice(1:25) %>%
  # filter(start < end) %>%
  ungroup -> sel

# Load corresponding genomes
sel %>%
  mutate(path = sprintf('data/A_representatives/%s/genome.fna.gz',
                        tax.bio)) %>%
  pull(path) %>%
  unique %>%
  future_map(readDNAStringSet) %>%
  invoke(.f = c) -> sel.genomes

# Extract sequence as determined by the computations above
sel %>%
  select(seqnames = tax.bio.chr, start, end, strand) %>%
  as_granges() -> sel.range

sel$observed <- BSgenome::getSeq(sel.genomes, sel.range) %>% as.character()

# Load sequences listed in files
seq.helper <- function(path) {
  path %>%
    # read in Stockholm
    read_lines %>%
    # Keep lines with nucleotide info
    keep(str_detect, '^[0-9]') %>%
    tibble(row = .) %>%
    # separate gene id and sequence
    separate(row, c('gene', 'seq'), sep = ' +') %>%
    mutate_at('seq', str_replace_all, 'U', 'T') %>%
    mutate_at('seq', str_remove_all, '[^A-Z]') %>%
    # keep track of position of lines
    mutate(r = 1:dplyr::n()) %>%
    # combine for genes with info in multiple lines
    # in order of apperane in file
    arrange(gene, r) %>%
    group_by(gene) %>%
    summarize(stockholm.seq = str_c(seq, collapse = '')) %>%
    ungroup
}

# sel %>%
#   select(motif) %>%
#   unique %>%
#   left_join(motifs, 'motif') %>%
motifs %>%
  with(set_names(path, motif)) %>%
  future_map(seq.helper) %>%
  map2(names(.), ~ mutate(.x, motif = .y)) %>%
  bind_rows() -> expected


# Unforunately not given for now
# assertthat::are_equal(
#   expected %>%
#     right_join(sel, c('motif', 'gene' = 'tax.bio.gene')) %>%
#     filter(observed != stockholm.seq) %>%
#     nrow,
#   0
# )


################################################################################

# For now check for overlaps with Rfam matches for those that were directly 
# lookup-able

expected %>%
  right_join(sel, c('motif', 'gene' = 'tax.bio.gene')) %>%
  filter(observed == stockholm.seq) %>%
  select(
    motif, region, side, tax.bio,
    seqnames = tax.bio.chr,
    start, end, strand
  ) -> xs



xs %>%
  select(motif, tax.bio) %>%
  unique %>%
  count(tax.bio) %>%
  ggplot(aes(n)) +
  geom_histogram() +
  xlab('Motifs per tax.bio')

xs %>%
  count(motif, tax.bio) %>%
  ggplot(aes(n)) +
  geom_histogram() +
  scale_x_log10() +
  xlab('No of motif copies in species')



xs %>%
  select(motif, tax.bio) %>%
  unique %>%
  count(motif) %>%
  ggplot(aes(n)) +
  geom_histogram() +
  xlab('Different species in motif')
  

################################################################################

rfam <- read_tsv(in.rfam)

rfam %>%
  mutate(
    tmp = start,
    start = ifelse(strand == '+', start, end),
    end = ifelse(strand == '+', end, tmp)
  ) %>%
  select(seqnames = chr, start, end, strand,
         rf, name) %>%
  plyranges::as_granges() %>%
  plyranges::mutate(rf.len = width, rf.strand = strand) -> rf.range


xs %>%
  plyranges::as_granges() %>%
  plyranges::mutate(motif.len = width, motif.strand = strand) -> motif.range



inter <- plyranges::join_overlap_intersect(motif.range, rf.range)


inter %>%
  as_tibble %>%
  count(rf.strand == motif.strand)

inter %>%
  as_tibble %>%
  mutate(shared = width) %>%
  ggplot(aes(shared /rf.len, shared / motif.len)) +
  geom_point() -> p

ggExtra::ggMarginal(p)

inter %>%
  as_tibble() %>%
  filter(width / motif.len > .95) %>%
  select(motif) %>%
  unique -> recalling
  # count(rf, name)

################################################################################
# swoop in for SYP PCC 7002

expected %>%
  filter(str_starts(gene, '32049.SAMN01081740')) %>%
  select(motif) %>%
  unique -> recalling
################################################################################

stats <- read_tsv('data/I_fdr.tsv')

stats %>%
  left_join(mutate(recalling, foo = TRUE)) %>%
  mutate_at('foo', replace_na, FALSE) %>%
  select(- c(region, region.avgid, region.gc)) %>%
  filter(foo) %>%
  filter(
    RNAphylo.fdr < 10,
    alignment.power > 10,
    fraction.covarying > 10,
    hmmpair > 1
  ) %>%
  View
  
  head
  View
  gather('k', 'value', -c(motif, foo)) %>%
  dplyr::rename('overlap Rfam > 95% at least once' = foo) %>%
  # filter(k == 'RNAphylo.fdr') %>%
  # mutate_at('value', log10) %>%
  ggplot(aes(value, color = `overlap Rfam > 95% at least once`)) +
  stat_ecdf() +
  facet_wrap(~ k , scales = 'free')
    

################################################################################
