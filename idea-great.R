library(tidyverse)
library(rGREAT)

library(doParallel)
library(furrr)

cores <- 8
plan(multisession, workers = cores)

################################################################################
# Input data from project

'data/J_novel/all_intergenic_regions.tsv.gz' |>
  read_tsv() -> regions

'data/J_novel/references_inside-of_intergenic_regions.tsv.gz' |>
  read_tsv() -> ref.regions

################################################################################
# Get genome lengths

# Cannot directly use
# 'data/A_checkm/checkm_summary.tsv'
# because summarized for entire organism, not by chromosome

# Thus: Manually check fasta files
my_genome_lens <-
  'data/A_representatives/*/genome.fna.gz' %>%
  Sys.glob() |>
  future_map(Biostrings::fasta.seqlengths) |>
  unlist()


################################################################################
# Load curated KEGG annotation

kegg.path.ko <- 'https://rest.kegg.jp/link/pathway/ko' |>
  read_tsv(col_names = c('term', 'path')) |>
  mutate_all(str_remove, '^.*:')

kegg.paths <- 'https://rest.kegg.jp/list/pathway' |>
  read_tsv(col_names = c('path', 'pathway'))

kegg.ko <- 'https://rest.kegg.jp/list/ko' |>
  read_tsv(col_names = c('term', 'ortholog'))

################################################################################
# Prepare for genomic range enrichment

# Input, regions of interest from Rfam hits
# As list per Rfam entry
gr.input <-
  ref.regions |>
  filter(type == 'Rfam') |>
  mutate(
    name = sprintf('%s (%s)', name, rfam)
  ) |>
  select(seqnames, start, end, strand, name) |>
  group_by(name) |>
  do(gr = . |>
       select(- name) |>
       plyranges::as_granges()
       ) |>
  with(set_names(gr, name))

# Background, all search regions
regions2 <-
  regions |>
  mutate(row = paste0('region-row-', 1:n())) 

gr.bg <-
  regions2 |>
  select(
    seqnames, start, end, strand,
    row
  ) |>
  plyranges::as_granges()

# Foreground "TSS", aka search regions with hits
gr.tss <-
  regions2 |>
  semi_join(ref.regions, 'region') |>
  mutate(
    tss_strand = strand,
    strand = '*',
    tss_position = map2(start, end, mean) |> unlist()
  ) |>
  select(
    seqnames, start, end, strand,
    gene_id = row, tss_position, tss_strand
  ) |>
  plyranges::as_granges()

# Custom region-pathway geneset map
gr.gs <-
  regions2 |>
  select(region, row) |>
  separate(region, c('term', 'side'), sep = '_') |>
  left_join(kegg.path.ko, 'term', relationship = "many-to-many") |>
  inner_join(kegg.paths, 'path') |>
  unite('p2', path, pathway, sep = ': ') |>
  select(p2, row) |>
  unique() |>
  group_by(p2) |>
  do(i = list(.$row)) |>
  ungroup() |>
  with(set_names(i, p2)) |>
  map(1)

################################################################################
# Custom GREAT run

source('scripts/helper-great.R')

great.res <-
  gr.input |>
  # rather parallize per family
  futrue_map(function(i) {
    my.great(
      gr = i,
      gene_sets = gr.gs,
      my_genome_lens = my_genome_lens,
      extended_tss = gr.tss,
      background = gr.bg,
      cores = 1
    )
  })


################################################################################
getEnrichmentTable(foo) |>
  filter(p_adjust_hyper <= 0.05)

jpeg('foo.jpg')
rGREAT::plotVolcano(foo)
dev.off()


