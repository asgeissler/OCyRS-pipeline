# Purpose:
# Look up the coordinates of the detect motifs in the bacterial genomes

library(tidyverse)
library(Biostrings)
library(plyranges)
library(furrr)

in.path <- 'data/I_fdr.tsv'
in.genes <- 'data/A_representatives/genes.tsv.gz'

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

# Nuisance: Determining positions for 'left' cases (see pos.helper comments)
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
#
# Case:      Left                Right         Left                    Right

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
      (side == 'upstream')   & (gene.strand == '+') ~ gene.start - region.len + rel.start,
      (side == 'upstream')   & (gene.strand == '-') ~ gene.end   + region.len - rel.end,
      (side == 'downstream') & (gene.strand == '+') ~ gene.end.               + rel.start,
      (side == 'downstream') & (gene.strand == '-') ~ gene.start - region.len + rel.end
    ),
    end = case_when(
      (side == 'upstream')   & (gene.strand == '+') ~ gene.start - region.len + rel.end,
      (side == 'upstream')   & (gene.strand == '-') ~ gene.end   + region.len - rel.start,
      (side == 'downstream') & (gene.strand == '+') ~ gene.end.               + rel.end,
      (side == 'downstream') & (gene.strand == '-') ~ gene.start - region.len + rel.start
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
    select(- starts_with('rel'), - starts_with('gene'), - is.left)  %>%
    mutate(motif = motif)
}


motifs %>%
  select(path, side, region, motif) %>%
  head %>%
  future_pmap(worker) %>%
  bind_rows() -> abs.pos


################################################################################
################################################################################
# Though not strictly needed, select some upstream/downstream cases
# and double check that the computed position are matching the sequences.
# Doing such a test ensures that the assumption on CMfinder's output is given.

# select 5 arbitrary examples for each scenario (see comments pos.helper)
abs.pos %>%
  mutate(side = ifelse(str_detect(motif, '_upstream.fna.motif'),
                       'upstream', 'downstream')) %>%
  group_by(side, strand) %>%
  slice(1:5) %>%
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

sel$observed <- getSeq(sel.genomes, sel.range) %>% as.character()

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

sel %>%
  select(motif) %>%
  unique %>%
  left_join(motifs, 'motif') %>%
  group_by(motif, path) %>%
  do(i = seq.helper(.$path)) %>%
  ungroup %>%
  unnest(i) -> expected

expected %>%
  right_join(sel, c('motif', 'gene' = 'tax.bio.gene')) %>%
  # select(gene, motif, observed, stockholm.seq) %>%
  filter(observed != stockholm.seq) %>%
  View
sel
expected
################################################################################