# Code based on `great` function from rGREAT
# https://github.com/jokergoo/rGREAT/blob/master/R/great_local.R
# Purpose: Code adapted to avoid dependency on GenomeInfoDb 
# Needed for support of non human/mouse/etc genomes

my.great <- function(
    gr, gene_sets, extended_tss, background,
    my_genome_lens, # should be named integer vector
    min_gene_set_size = 5,
    cores = 1
) {
  param = list()
  gene_sets_name = ""
  
  mt = attributes(metadata(extended_tss))
  param$genome = "unknown"
  param$tss_source = "self-provided"
  param$orgdb = 'NA'
  
  extended_tss = sort(extended_tss)
  mode_param = attr(extended_tss, "mode")
  param = c(param, mode_param)
  
  mt = attributes(metadata(extended_tss))
  gene_sets_name = "self-provided"
  
  gene_sets = lapply(gene_sets, unique)
  names(extended_tss) = extended_tss$gene_id
  
  message("* check gene ID type in `gene_sets` and in `extended_tss`.")
  
  v1 = unique(unlist(gene_sets))
  if(length(intersect(v1, extended_tss$gene_id))/length(v1) <= 0) {
    stop_wrap(qq("It seems the gene ID type in `gene_sets` (e.g. '@{gene_sets[[1]][1]}') is different from in `extended_tss` (e.g. '@{names(extended_tss)[1]}')."))
  }
  
  sl = my_genome_lens
  gr_genome = GRanges(seqnames = names(sl), ranges = IRanges(1, sl))
  background_type = "whole genome"
  if(is.null(background)) {
    background = gr_genome
    message("* use whole genome as background.")
  } else {
    message("* use self-defined background regions.")
    if(inherits(background, "character")) {
      background_chr = background
      if(!is.null(biomart_dataset)) {
        background_chr = gsub("^chr", "", background_chr)
      }
      background = gr_genome[seqnames(gr_genome) %in% background_chr]
    }
    background_type = "self-provided"
  }
  strand(background) = "*"
  
  # skip exclude part
  
  background = reduce(background)
  background_total_length = sum(width(background))
  
  strand(gr) = "*"
  
  gr_origin = gr
  
  gr_mid = gr
  start(gr_mid) = end(gr_mid) = mid(gr)
  gr = gr_mid
  
  message("* overlap `gr` to background regions (based on midpoint).")
  
  gr = intersect(gr, background)  # possible seqlevels are different due to different versions of organisms
  
  n_total = length(gr)
  
  message(GetoptLong::qq("* in total @{n_total} `gr`."))
  message(GetoptLong::qq("* overlap extended TSS to background regions."))
  
  ov2 = findOverlaps(extended_tss, background)
  r1 = extended_tss[queryHits(ov2)]
  r2 = background[subjectHits(ov2)]
  
  # note in extended_tss, there might be same gene in multi rows, due to a gene is segmented by background
  extended_tss2 = pintersect(r1, r2)
  extended_tss2$gene_id = r1$gene_id
  names(extended_tss2) = r1$gene_id
  extended_tss = extended_tss2
  extended_tss$hit = NULL
  
  # note, after overlap to background, an extended TSS may be split into several regions
  message("* check which genes are in the gene sets.")
  gene_sets = lapply(gene_sets, function(x) intersect(x, names(extended_tss)))
  
  message(GetoptLong::qq("* only take gene sets with size >= @{min_gene_set_size}."))
  l = sapply(gene_sets, length) >= min_gene_set_size
  gene_sets = gene_sets[l]
  gene_sets_ind = lapply(gene_sets, function(x) which(names(extended_tss) %in% x))
  
  n_set = length(gene_sets)
  
  if(n_set == 0) {
    stop_wrap("No gene set left.")
  }
  
  message(GetoptLong::qq("* in total @{n_set} gene sets."))
  
  p = numeric(n_set)
  n_obs = numeric(n_set)
  n_exp = numeric(n_set)
  fold = numeric(n_set)
  prop = numeric(n_set)
  mean_tss_dist = numeric(n_set)
  
  gene_hits = numeric(n_set)
  
  message(GetoptLong::qq("* overlap `gr` to every extended TSS."))
  ov = findOverlaps(gr, extended_tss)
  
  message("* perform binomial test for each biological term.")
  pb = progress::progress_bar$new(total = n_set)
  
  all_chr = as.vector(seqnames(extended_tss))
  all_s = start(extended_tss)
  all_e = end(extended_tss)
  tss_position = extended_tss$tss_position
  gr_start = start(gr)
  
  registerDoParallel(cores)
  lt <- foreach(i = seq_len(n_set)) %dopar% {
    if(cores <= 1) {
      pb$tick()
    }
    ind = gene_sets_ind[[i]]
    
    if(length(ind) == 0) {
      width_fgr = 0
    } else {
      sl = split(all_s[ind], all_chr[ind])
      el = split(all_e[ind], all_chr[ind])
      width_fgr = sum(sapply(seq_along(sl), function(i) {
        reduce_by_start_and_end(sl[[i]], el[[i]])
      }))
    }
    prop = width_fgr/background_total_length
    
    l = subjectHits(ov) %in% ind
    n_hits = length(unique(queryHits(ov)[l]))
    
    mean_tss_dist = mean(abs(tss_position[subjectHits(ov)[l]] - gr_start[queryHits(ov)[l]]))
    
    if(n_hits == 0) {
      p = 1
    } else {
      p = 1 - pbinom(n_hits - 1, n_total, prop)
    }
    n_obs = n_hits
    n_exp = prop*n_total
    
    gene_hits = length(unique(extended_tss$gene_id[subjectHits(ov)[l]]))
    
    return(list(prop = prop, n_obs = n_obs, n_exp = n_exp, p = p, mean_tss_dist = mean_tss_dist, gene_hits = gene_hits))
  }
  stopImplicitCluster()
  
  prop = sapply(lt, function(x) x$prop)
  n_obs = sapply(lt, function(x) x$n_obs)
  n_exp = sapply(lt, function(x) x$n_exp)
  p = sapply(lt, function(x) x$p)
  mean_tss_dist = sapply(lt, function(x) x$mean_tss_dist)
  gene_hits = sapply(lt, function(x) x$gene_hits)
  
  fold = n_obs/n_exp
  
  df = data.frame(id = names(gene_sets),
                  genome_fraction = prop,
                  observed_region_hits = n_obs,
                  fold_enrichment = fold,
                  p_value = p,
                  p_adjust = p.adjust(p, "BH"),
                  mean_tss_dist = round(mean_tss_dist),
                  observed_gene_hits = gene_hits,
                  gene_set_size = sapply(gene_sets, length))
  
  n_gene_total = length(unique(extended_tss$gene_id))
  n_gene_gr = length(unique(extended_tss[subjectHits(ov)]$gene_id))
  df$fold_enrichment_hyper = df$observed_gene_hits/(df$gene_set_size*n_gene_gr/n_gene_total)
  df$p_value_hyper = 1 - phyper(df$observed_gene_hits - 1, df$gene_set_size, n_gene_total - df$gene_set_size, n_gene_gr)
  df$p_adjust_hyper = p.adjust(df$p_value_hyper, "BH")
  
  df = df[order(df$p_adjust, df$p_value, -df$fold_enrichment), , drop = FALSE]
  
  if(all(grepl("^GO:\\d+$", df$id))) {
    check_pkg("GO.db", bioc = TRUE)
    terms = Term(GO.db::GOTERM)
    
    desc = unname(terms[df$id])
    desc[is.na(desc)] = ""
    df = cbind(data.frame(id = df$id, description = desc),
               df[, -1])
  }
  rownames(df) = NULL
  
  obj = GreatObject(
    table = df,
    gr = gr_origin,
    n_total = n_total,
    gene_sets = gene_sets,
    gene_sets_name = gene_sets_name,
    extended_tss = extended_tss,  # has been intersected with background,
    n_gene_gr = n_gene_gr,
    background = background,
    background_type = background_type, 
    param = param
  )
  
  return(obj)
}
