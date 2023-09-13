# Prepare output for PETfold
#
# Part 1:
# Adapted from J_aln-seq-pos.R
# Convert Stockholm alignment to FastA (including gaps)
# Part 2:
# - export secondary structure
# - for motifs, also re-create tree for only used labels

library(tidyverse)
library(ape)
library(furrr)

# Assume exists, the motifs
#  data/J_novel/potentially-novel-motifs.tsv
# Rfam files, eg with implicit dependency via
# data/G_rfam-cmsearch.tsv.gz

cpus <- as.integer(unlist(snakemake@threads))
# cpus <- 20
plan(multisession, workers = cpus)


###############################################################################
# Helper for reading/writing

worker <- function(in.sto,  out.aln, out.struct, in.tree = NULL, out.tree = NULL) {
  # in.sto <- 'data/G_rfam-bacteria-seeds/RF00001.sto'
  # in.sto <- 'data/F_cmfinder/D_search-seqs/K00003_upstream/K00003_upstream.fna.motif.h1_1'
  # in.tree <- 'data/C_shrink/K00003/output.tree'
  
  dir.create(dirname(out.aln))
  
  # read in Stockholm
  dat <- read_lines(in.sto)
  
  # Process alignment
  aln <-
    dat |>
    # Keep lines with nucleotide info
    keep(str_detect, '^[A-Za-z0-9]') %>%
    tibble(row = .) %>%
    # separate gene id and sequence
    separate(row, c('gene', 'seq'), sep = ' +') %>%
    # !Difference to J_aln-seq-pos.R -> Keep gaps etc
    # mutate_at('seq', str_replace_all, 'U', 'T') %>%
    # mutate_at('seq', str_remove_all, '[^A-Z]') %>%
    # keep track of position of lines
    mutate(r = 1:n()) %>%
    # combine for genes with info in multiple lines
    # in order of appearance in file
    arrange(gene, r) %>%
    group_by(gene) %>%
    summarize(stockholm.seq = str_c(seq, collapse = '')) %>%
    ungroup %>%
    with(set_names(stockholm.seq, gene)) |>
    Biostrings::RNAStringSet()
  
  Biostrings::writeXStringSet(aln, out.aln)
  
  # Load secondary structure
  dat |>
    keep(str_starts, '#=GC SS_cons') |>
    str_remove('^#=GC SS_cons *') |>
    str_c(collapse = '') |>
    str_replace_all('<', '(') |>
    str_replace_all('>', ')')  |>
    str_replace_all('[^()]', '.') |>
    write_lines(out.struct)
  
  # Process tree
  if (!is.null(in.tree)) {
    x <-
      in.tree |>
      read.tree() |>
      # Only keep tips with sequence in alignment
      keep.tip(names(aln))
    
    write.tree(x, out.tree)
  }
}

###############################################################################
# 1. Work on all bacterial Rfam families

dir.create('data/M_alignments-Rfam/')

tibble(
  in.sto = Sys.glob('data/G_rfam-bacteria-seeds/*.sto'),
  rf = in.sto |>
    basename() |>
    str_remove('\\.sto'),
  out.aln =    sprintf('data/M_alignments-Rfam/%s/aln.fna', rf),
  out.struct = sprintf('data/M_alignments-Rfam/%s/structure.txt', rf)
) |>
  select(- rf) |>
  future_pwalk(worker)

###############################################################################
# 2. Work on the ~400 candidate novel motifs

dir.create('data/M_alignments-candidates/')

tibble(
  motif = 'data/J_novel/potentially-novel-motifs.tsv' |>
    read_tsv() |>
    pull(motif),
  region = str_remove(motif, '.fna.motif.*'),
  ortho  = str_remove(region, '_[updownstream]+$'),
  in.sto     = sprintf('data/F_cmfinder/D_search-seqs/%s/%s', region, motif),
  in.tree    = sprintf('data/C_shrink/%s/output.tree', ortho),
  out.aln    = sprintf('data/M_alignments-candidates/%s/aln.fna', motif),
  out.struct = sprintf('data/M_alignments-candidates/%s/structure.txt', motif),
  out.tree   = sprintf('data/M_alignments-candidates/%s/tree.txt', motif)
) |>
  select(- c(motif, region, ortho)) |>
  future_pwalk(worker)

