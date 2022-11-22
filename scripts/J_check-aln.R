# Load sequences of motif hits
# Check if these are part of the motif alignments
# Determine per motif score cutoff

library(tidyverse)
library(furrr)


# in.cmsearch <- 'data/J_cmsearch-collected.tsv'
# in.motifs <- 'data/F_cmfinder/motifs-demerged.tsv'
in.cmsearch <- unlist(snakemake@input[['collected']])
in.motifs <- unlist(snakemake@input[['motifs']])

# out.motifpos <- 'data/J_motif-positions.tsv'
# out.cutoffs <- 'data/J_score-cutoffs.tsv'
# out.homologs <- 'data/J_motif-homologs.tsv'
# out.fig <- 'data/J_overlaps.png'

out.motifpos <- unlist(snakemake@output[['motifpos']])
out.cutoffs <- unlist(snakemake@output[['cutoffs']])
out.homologs <- unlist(snakemake@output[['homologs']])
out.fig <- unlist(snakemake@output[['fig']])


cpus <- as.integer(unlist(snakemake@threads))
# cpus <- 20
plan(multisession, workers = cpus) 

###############################################################################
# Load coordinates of hits

cms <- read_tsv(in.cmsearch)

# Split per genome
cms %>%
  rename(motif = name) %>%
  mutate(
    tmp = start,
    start = ifelse(start < end, start, end),
    end = ifelse(tmp < end, end, tmp)
  ) %>%
  select(- tmp) %>%
  group_by(tax_bio) %>%
  do(i = {
    grp <- .
    grp %>%
      select(seqnames = tax.bio.chr, start, end, strand,
             motif, score, evalue) %>%
      plyranges::as_granges()
  }) %>%
  ungroup %>%
  with(set_names(i, tax_bio)) -> cms.ranges


###############################################################################
# Load sequences of alignments

worker <- function(path, motif) {
  # path <- 'data/F_cmfinder/D_search-seqs/K00003_upstream/K00003_upstream.fna.motif.h1_1'
  # Load the relative positions from the Stockholm motif file
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
    mutate(r = 1:n()) %>%
    # combine for genes with info in multiple lines
    # in order of apperance in file
    arrange(gene, r) %>%
    group_by(gene) %>%
    summarize(stockholm.seq = str_c(seq, collapse = '')) %>%
    ungroup %>%
    mutate(tax.bio = str_remove(gene, '\\.[^.]*$')) %>%
    select(tax.bio, stockholm.seq) %>%
    mutate(motif = motif)
}

in.motifs %>%
  read_tsv() %>%
  filter(dir == 'D_search-seqs') %>%
  select(path, filename) %>%
  with(future_map2(path, filename, worker)) %>% 
  bind_rows() -> motif.seqs
  

###############################################################################
# Find locations of all aligned sequences following the
# Efficient genome searching tutorial, section
# Finding an arbitrary nucleotide pattern in an entire genome
# https://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/GenomeSearching.pdf
# Idea: reverse complement pattern before search

motif.seqs$tax.bio %>%
  unique %>%
  future_map(function(i) {
    sprintf('data/A_representatives/%s/genome.fna.gz', i) %>%
      Biostrings::readDNAStringSet() -> genome
    
    motif.seqs %>%
      filter(tax.bio == i) %>%
      mutate(row = 1:n()) %>%
      mutate(strand = '+') -> fwd
    
    fwd %>%
      mutate(
        stockholm.seq = stockholm.seq %>%
          Biostrings::DNAStringSet() %>%
          Biostrings::reverseComplement() %>%
          as.character(),
        strand = '-'
      ) %>%
      bind_rows(fwd) -> patterns
    
    
    patterns %>%
      group_by_all() %>%
      do(search = {
        grp <- .
        Biostrings::vmatchPattern(grp$stockholm.seq, genome) %>%
          as.data.frame() %>%
          transmute(start, end, chr = names(genome)[group])
      }) %>%
      ungroup %>%
      unnest(search) %>%
      select(motif, chr, start, end, strand)
  }) %>%
  bind_rows() -> motif.seqs.pos

write_tsv(motif.seqs.pos, out.motifpos)

###############################################################################
# Check overlap between hits and alignment sequences

names(cms.ranges) %>%
  future_map(function(i) {
    motif.seqs.pos %>%
      filter(str_starts(chr, i)) %>%
      rename(seqnames = chr) %>%
      plyranges::as_granges() %>%
      mutate(aln.len = IRanges::width(.)) -> motif.range
    
    cms.ranges[[i]] %>%
      mutate(hit.len = IRanges::width(.)) %>%
      plyranges::join_overlap_intersect_directed(motif.range) %>%
      as_tibble() %>%
      filter(motif.x == motif.y)
  }) %>%
  bind_rows() %>%
  mutate(jacc = width / (hit.len + aln.len - width)) -> cmsearch.motif.overlap

###############################################################################

cmsearch.motif.overlap %>%
  ggplot(aes(jacc)) +
  stat_ecdf() +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  xlab('Jaccard similarity cmsearch hit\nand motif alignment seq. pos.') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  ylab('Cum. emp. density') +
  theme_bw(18) -> p1

###############################################################################
cmsearch.motif.overlap %>%
  filter(jacc > .99) %>%
  group_by(motif = motif.x) %>%
  summarize(min.seq.score = min(score)) %>%
  ungroup -> cutoff

write_tsv(cutoff, out.cutoffs)

###############################################################################

cms %>%
  left_join(cutoff, c('name' = 'motif')) %>%
  drop_na -> cand

cand %>%
  ggplot(aes(score / min.seq.score, - log10(evalue))) +
  geom_hex() +
  scale_fill_viridis_c() +
  xlab('ratio cmsearch hit score / \n max. motif alignment score') +
  theme_bw(18) -> p2

cmsearch.motif.overlap %>%
  filter(jacc > .99) %>%
  # pull(evalue) %>% summary
  ggplot(aes(evalue)) +
  geom_histogram(bins = 50) +
  scale_x_log10(breaks = c(1e-40, 1e-20, 1e-10, .05)) +
  xlab('max. E-value alignment score per motif') +
  theme_bw(18) -> p3

###############################################################################

cand %>%
  filter(score > min.seq.score) %>%
  filter(evalue < .05) -> homologs

write_tsv(homologs, out.homologs)

###############################################################################

full_join(
  homologs %>%
    count(motif = name, tax_bio) %>%
    group_by(motif) %>%
    summarize(avg.homo = mean(n)),
  motif.seqs %>%
    count(motif) %>%
    rename(no.seq = n),
  'motif'
) %>%
  ggplot(aes(avg.homo)) +
  stat_ecdf() +
  scale_x_log10() +
  xlab('Average no. predicted homolog motifs per genome\nafterwith E-value 0.05 and motif score') +
  annotation_logticks(sides = 'b') +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  ylab('Cum. emp. density') +
  theme_bw(18) -> p4

###############################################################################

cowplot::plot_grid(
  p1, p2, p3, p4,
  labels = 'AUTO',
  label_size = 18
)

ggsave(out.fig, width = 20, height = 12)

###############################################################################