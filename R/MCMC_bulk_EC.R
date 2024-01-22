MCMC_bulk_EC = function(SE,
                        EC_data,
                        min_counts_ECs,
                        n_samples,
                        n_samples_per_group,
                        numeric_groups,
                        sample_ids_per_group,
                        n_groups,
                        tr_ids_sce,
                        tr_ids_sce_pass_filter,
                        levels_groups,
                        N_MCMC,
                        burn_in,
                        n_cores,
                        undersampling_int,
                        traceplot){
  c_prop = 0.3 # proportionality constant of the ARW
  
  counts = EC_data[[1]]
  list_EC_tr_id = EC_data[[2]]
  list_EC_US_id = EC_data[[3]]
  tr_id = EC_data[[4]]
  
  tot_n_counts = sum(unlist(counts))
  
  n_tr = length(tr_id)
  
  # ok!
  #n_tr;
  #min(unlist(list_EC_tr_id)); max(unlist(list_EC_tr_id))
  
  rm(EC_data)
  
  S = assays(SE)$spliced
  U = assays(SE)$unspliced
  
  has_U_tr = rowData(SE)$has_U_tr
  
  eff_len_S = rowData(SE)$eff_len_S
  eff_len_U = rowData(SE)$eff_len_U
  
  #  max_len = max(c(eff_len_S, eff_len_U), na.rm = TRUE)
  #  eff_len_S = eff_len_S/max_len
  #  eff_len_U = eff_len_U/max_len
  
  rm(SE)
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # separate uniquely mapping ECs:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  list_X_unique = list()
  
  # separate uniquely mapping ECs:
  for(sample in seq_len(n_samples)){
    # fill X with 0's (later replace with uniquely mapping reads):
    X = matrix( 0, nrow = n_tr, ncol = 2)
    
    # unique bool vector indicating uniquely mapping ECs.
    n_EC = length(list_EC_tr_id[[sample]])
    unique = rep(FALSE, n_EC)
    
    SU_id = sel = c()
    
    for(ec in seq_len(n_EC)){
      gene_id = list_EC_tr_id[[sample]][[ec]] + 1
      if(length(gene_id) == 1){
        unique[ec] = TRUE
        SU_id = list_EC_US_id[[sample]][[ec]] + 1
        sel = cbind(gene_id, SU_id)
        
        X[ sel ] = X[ sel ] + counts[[sample]][ec]
      }
    }
    
    list_X_unique[[sample]] = X
    
    # remove "unique" ECs from latent variables (already assigned in "list_X_unique")
    list_EC_tr_id[[sample]] = list_EC_tr_id[[sample]][!unique]
    list_EC_US_id[[sample]] = list_EC_US_id[[sample]][!unique]
    counts[[sample]] = counts[[sample]][!unique]
  }
  
  #n_tr;
  #min(unlist(list_EC_tr_id)); max(unlist(list_EC_tr_id))
  #dim(list_X_unique[[1]])
  
  rm(X); rm(unique); rm(sel); rm(SU_id); rm(n_EC)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # FILTER 1) remove ECs with 0 counts (within a specific cell type):
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(min_counts_ECs > 0){
    
    # 1) remove ECs with 0 counts (within a specific cell type):
    sel_non_zero_EC = lapply(counts, function(X){
      X > min_counts_ECs
    })
    
    # reduction in ECs:
    xx = 1-sum(sapply(sel_non_zero_EC, sum))/sum(sapply(sel_non_zero_EC, length))
    
    # reduction in counts:
    yy = 1-sum(sapply(seq_len(n_samples), function(id){
      sum(counts[[id]][sel_non_zero_EC[[id]]])}))/tot_n_counts
    
    message(xx*100, " % of ECs discarded due to 'min_counts_ECs' filter.")
    message(yy*100, " % of counts discarded  due to 'min_counts_ECs' filter.")
    
    # trim EC gene list, each element = 1 EC, OK:
    list_EC_tr_id = lapply(seq_len(n_samples), function(X){
      list_EC_tr_id[[X]][ sel_non_zero_EC[[X]] ]
    })
    
    # trim EC SU list, each element = 1 EC, OK:
    list_EC_US_id = lapply(seq_len(n_samples), function(X){
      list_EC_US_id[[X]][ sel_non_zero_EC[[X]] ]
    })
    
    # trim counts, each element = 1 EC, OK:
    counts = lapply(seq_len(n_samples), function(X){
      counts[[X]][ sel_non_zero_EC[[X]] ]
    })
    
    rm(sel_non_zero_EC)
  }
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # FILTER 2) remove tr_id absent across ALL ECs of ALL tr_id...check SE of those tr_id -> it'd be ~ 0!
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # Select tr_id in ECs with non-zero counts:
  tr_non_zero = tr_id[ unique(unlist(list_EC_tr_id)) + 1 ]
  # PLUS ADD GENES with UNIQUELY MAPPING READS!!!
  tr_non_zero = c(tr_non_zero,  tr_id[ rowSums(do.call(cbind, list_X_unique)) > 0 ] )
  tr_non_zero = unique(tr_non_zero)
  tr_non_zero = sort(tr_non_zero)
  
  # tr_id
  matches = match(tr_id, tr_non_zero)
  
  list_EC_tr_id = lapply(list_EC_tr_id, function(X){
    lapply(X, function(x){
      matches[x + 1] - 1
    })
  })
  
  matches = match(tr_non_zero, tr_id)
  list_X_unique = lapply(list_X_unique, function(x){
    # need to do: match(tr_id, tr_non_zero)
    # tr_non_zero don't necessarily follow the same order of tr_id!!
    x[matches,] #
  })
  
  rm(matches)
  
  n_tr_non_zero = length(tr_non_zero)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # set starting values, defined from "SE" values
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  TOT_counts = S + U
  # tr_ids_sce is the row name of TOT_counts
  
  # use pi_gene, estimated from sample-specific counts (+ 1):
  matches = match(tr_non_zero, tr_ids_sce)
  
  # match effective lengths:
  eff_len_S = eff_len_S[matches]
  eff_len_U = eff_len_U[matches]
  # if no match -> set to 1? maybe to the average length?
  nas = is.na(eff_len_S)
  eff_len_S[nas] = eff_len_U[nas] = mean(eff_len_S, na.rm = TRUE)
  # if NA's still present -> replace with 1
  nas = is.na(eff_len_S)
  eff_len_S[nas] = eff_len_U[nas] = 1
  
  PI_gene = 1 + TOT_counts[matches, ]
  rownames(PI_gene) = tr_non_zero
  PI_gene[is.na(PI_gene)] = 1 # set NAs to 1
  PI_gene[ PI_gene < 1 ] = 1 # set counts < 1 to 1
  PI_gene = apply(PI_gene, 2, function(x) x/sum(x))
  
  # assign pi_SU as starting MCMC value:
  # for transcripts without U version (has_U_tr == FALSE), set PI_SU to (1,0)
  
  tot_spliced = 1/2 + S
  # tr without U version, set S to 1 and U to 0:
  tot_spliced[ !has_U_tr, ] = 1
  tot_spliced = tot_spliced[matches, ]
  tot_spliced[is.na(tot_spliced)] = 1/2
  
  tot_unspliced = 1/2 + U
  # tr without U version, set S to 1 and U to 0:
  tot_unspliced[ !has_U_tr, ] = 0
  tot_unspliced = tot_unspliced[matches, ]
  tot_unspliced[is.na(tot_unspliced)] = 1/2
  
  PI_SU = lapply(seq_len(n_samples), function(x){
    X = cbind(tot_spliced[,x], tot_unspliced[,x])
    X/rowSums(X)
  })
  
  rm(matches)
  
  # transcripts to analyze: those that:
  # 1) appear in at least 1 non-zero EC (tr_non_zero)
  # 2) were not filtered in `SE` (min counts respected + have U version): tr_ids_sce_pass_filter
  # to guarantee that the order is preserved:
  tr_ids_keep = tr_non_zero[ tr_non_zero %in% tr_ids_sce_pass_filter ]
  
  keep_genes_id = which(tr_non_zero %in% tr_ids_keep)
  n_genes_keep = length(keep_genes_id)
  n_genes_keep
  
  rm(tr_non_zero)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # pi_S prior:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  keep_sce = tr_ids_sce %in% tr_ids_keep
  
  S = S[keep_sce,]
  U = U[keep_sce,]
  
  tr_ids_sce = tr_ids_sce[keep_sce]
  rm(keep_sce)
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # DRIMSeq prior:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  sel_non_zeros = (rowSums(S) > 0) & (rowSums(U) > 0)
  S = S[sel_non_zeros,]
  U = U[sel_non_zeros,]
  
  # compute global pi_S (across all samples):
  pi = rowSums(S)/( rowSums(S) + rowSums(U) )
  
  S_U = rbind(S, U)
  
  rm(S); rm(U)
  
  gene_2_tr = data.frame(gene_id = rep(tr_ids_sce[sel_non_zeros], 2),
                         transcript_id = c(paste(tr_ids_sce[sel_non_zeros], "S"), 
                                           paste(tr_ids_sce[sel_non_zeros], "U") ))
  rownames(S_U) = gene_2_tr$transcript_id
  
  # exp because `prior_precision` returns log-scale dispersion:
  log_precision = prior_precision(gene_to_transcript = gene_2_tr,
                                  transcript_counts = S_U,
                                  n_cores = n_cores,
                                  max_n_genes_used = 10^3)[[2]]
  
  rm(gene_2_tr); rm(S_U)
  
  prior_log_disp = c(mean(log_precision, na.rm = TRUE),
                     sd(log_precision, na.rm = TRUE))
  
  # PRIOR for log-delta_S:
  matches = match(tr_ids_sce[sel_non_zeros], names(log_precision))
  
  log_delta_S = log( pi * exp(log_precision[matches]))
  
  prior_log_S = c(mean(log_delta_S, na.rm = TRUE),
                  2*sd(log_delta_S, na.rm = TRUE))
  
  rm(log_delta_S); rm(log_precision); rm(pi);   rm(tr_ids_sce);  rm(sel_non_zeros); rm(matches)
  
  # if NA or NULL or Inf, use vaguely informative values:
  if( any(is.na(prior_log_disp) | is.infinite(prior_log_disp) | is.null(prior_log_disp)) ){
    prior_log_disp = c(log(2),5)
  }
  
  # if NA or NULL or Inf, use vaguely informative values:
  if( any(is.na(prior_log_S) | is.infinite(prior_log_S) | is.null(prior_log_S)) ){
    prior_log_S = c(log(2),5)
  }
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # initialize objects:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # TODO: if N_MCMC large -> use thinnings
  MCMC_bar_pi_1 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
  MCMC_bar_pi_2 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
  
  # delta_SU only for genes kept (n_genes_keep)
  delta_SU = lapply(seq_len(n_groups), function(g){
    ids = sample_ids_per_group[[g]] + 1
    n = length(ids)
    x = PI_SU[[ids[1]]][keep_genes_id,]
    if(n > 1){
      for(i in seq.int(2, n, by = 1)){
        x = x + PI_SU[[ids[i]]][keep_genes_id,]
      }
    }
    x = x/n * exp(prior_log_disp[1])
  })
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # MCMC:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  message("Starting the MCMC")
  
  PI_gene_times_SU = sapply(seq_len(n_samples), function(id){
    rep(PI_gene[,id], 2) * c(PI_SU[[id]][,1]/eff_len_S, PI_SU[[id]][,2]/eff_len_U )
  })
  #rm(PI_gene)
  
  list_X_unique = sapply(list_X_unique, c)
  
  # transform list_EC_tr_id, from 0 -> (n_tr-1) to 0 -> (2*n_tr-1)
  # if S: keep number in EC
  # if U: add n_tr_non_zero to number in EC
  list_EC_tr_id = lapply(seq_len(n_samples), function(sample){
    lapply(seq_along(list_EC_tr_id[[sample]]), function(ec){
      as.integer(list_EC_tr_id[[sample]][[ec]] + list_EC_US_id[[sample]][[ec]] * n_tr_non_zero)
    })
  }) # gene ids as in "gene_id" file.
  rm(list_EC_US_id)
  
  # only sample ECs every 10 iterations:
  sample_EC = rep(0, N_MCMC + 1)
  sample_EC[seq.int(1, N_MCMC, undersampling_int)] = 1
  
  sample_SU_TF = ifelse(seq_len(n_tr_non_zero) %in% keep_genes_id, 1, 0)
  
  res = .Call(`_DifferentialRegulation_Rcpp_MCMC_EC_US`,
              n_samples, # N samples
              n_tr_non_zero, # N tr_non_zero
              n_groups, # N groups
              n_genes_keep, # number of tr_non_zero to be analyzed.
              keep_genes_id-1, # vector indicating tr_non_zero to be analyzed (SU differential testing)
              numeric_groups - 1, # -1 ! # group id for every sample (must start from 0)
              sample_ids_per_group, # each list = vector with ids of samples 
              n_samples_per_group,
              N_MCMC, # MCMC iter
              burn_in, # burn-in
              PI_gene_times_SU,
              PI_SU, # prob of each gene (for every sample)
              list_X_unique, # SU uniquely mapping counts
              list_EC_tr_id, # SU uniquely mapping counts
              counts, # EC counts (integers)
              MCMC_bar_pi_1,
              MCMC_bar_pi_2,
              delta_SU,
              prior_log_disp,
              prior_log_S,
              sample_EC,
              sample_SU_TF,
              eff_len_S,
              eff_len_U,
              c_prop)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # check convergence:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  convergence = my_heidel_diag(res, R = N_MCMC, by. = 100, pvalue = 0.05)
  rm(res)
  
  # set convergence (over-written below if it did not converge)
  first = "converged"; second = "NOT run (1st chain converged)"
  
  if(convergence[1] == 0){ # if not converged, reset starting values and run a second chain (twice as long as the initial one):
    rm(convergence)
    
    message("Our MCMC did not converge (according to Heidelberger and Welch's convergence diagnostic):
                                       we will now double 'N_MCMC' and 'burn_in', and run it a second time.")
    
    N_MCMC = 2 * N_MCMC
    burn_in = 2 * burn_in
    
    # re-initialize objects:
    MCMC_bar_pi_1 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
    MCMC_bar_pi_2 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
    
    PI_SU = lapply(seq_len(n_samples), function(x){
      X = cbind(tot_spliced[,x], tot_unspliced[,x])
      X/rowSums(X)
    })
    
    PI_gene_times_SU = sapply(seq_len(n_samples), function(id){
      rep(PI_gene[,id], 2) * c(PI_SU[[id]][,1]/eff_len_S, PI_SU[[id]][,2]/eff_len_U )
    })
    
    rm(tot_spliced); rm(tot_unspliced); rm(PI_gene)
    
    delta_SU = lapply(seq_len(n_groups), function(g){
      ids = sample_ids_per_group[[g]] + 1
      n = length(ids)
      x = PI_SU[[ids[1]]][keep_genes_id,]
      if(n > 1){
        for(i in seq.int(2, n, by = 1)){
          x = x + PI_SU[[ids[i]]][keep_genes_id,]
        }
      }
      x = x/n * exp(prior_log_disp[1])
    })
    
    # only sample ECs every 10 iterations:
    sample_EC = rep(0, N_MCMC + 1)
    sample_EC[seq.int(1, N_MCMC, undersampling_int)] = 1
    
    res = .Call(`_DifferentialRegulation_Rcpp_MCMC_EC_US`,
                n_samples, # N samples
                n_tr_non_zero, # N tr_non_zero
                n_groups, # N groups
                n_genes_keep, # number of tr_non_zero to be analyzed.
                keep_genes_id-1, # vector indicating tr_non_zero to be analyzed (SU differential testing)
                numeric_groups - 1, # -1 ! # group id for every sample (must start from 0)
                sample_ids_per_group, # each list = vector with ids of samples 
                n_samples_per_group,
                N_MCMC, # MCMC iter
                burn_in, # burn-in
                PI_gene_times_SU,
                PI_SU, # prob of each gene (for every sample)
                list_X_unique, # SU uniquely mapping counts
                list_EC_tr_id, # SU uniquely mapping counts
                counts, # EC counts (integers)
                MCMC_bar_pi_1,
                MCMC_bar_pi_2,
                delta_SU,
                prior_log_disp,
                sample_EC,
                sample_SU_TF,
                eff_len_S,
                eff_len_U,
                c_prop)
    
    convergence = my_heidel_diag(res, R = N_MCMC, by. = 100, pvalue = 0.05)
    rm(res)
    
    if(convergence[1] == 0){ # if not converged for a 2nd time: return convergence error.
      first = "NOT converged"; second = "NOT converged"
      
      # create convergence DF:
      DF_convergence = data.frame(burn_in = NA,
                                  N_MCMC = N_MCMC,
                                  first_chain = first,
                                  second_chain = second)
      
      message("Our algorithm did not converged, try to increase N_MCMC.")
      
      return(list(NULL,
                  DF_convergence))
    }
    
    # if 2nd chain converged:
    first = "NOT converged"; second = "converged"
  }
  
  rm(list_EC_tr_id); #rm(list_EC_US_id)
  
  message("MCMC completed and successfully converged.")
  
  # the code below, is only run if either chain has converged:
  # increase the burn-in IF detected by "my_heidel_diag" (max burn-in = half of the chain length):
  burn_in = max(convergence[2]-1, burn_in)
  
  # create convergence DF:
  DF_convergence = data.frame(burn_in = burn_in,
                              N_MCMC = N_MCMC,
                              first_chain = first,
                              second_chain = second)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # compute p-value:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(traceplot){
    # store results of MCMC, before swapping:
    MCMC_U = MCMC_bar_pi_2
    names(MCMC_U) = levels_groups[1:2]
    MCMC_U$Transcript_id = tr_ids_keep
  }
  
  # discard burn-in
  sel = seq.int(from = burn_in +1, to = N_MCMC, by = 1)
  MCMC_bar_pi_1[[1]] = MCMC_bar_pi_1[[1]][sel,]
  MCMC_bar_pi_2[[1]] = MCMC_bar_pi_2[[1]][sel,]
  MCMC_bar_pi_1[[2]] = MCMC_bar_pi_1[[2]][sel,]
  MCMC_bar_pi_2[[2]] = MCMC_bar_pi_2[[2]][sel,]
  
  # swap A positions to decrease correlations:
  R = length(sel)
  swap = sample.int(R, R) # n indicates the nr of elements of the chain (exluded burn-in)
  MCMC_bar_pi_1[[1]] = MCMC_bar_pi_1[[1]][swap,]
  MCMC_bar_pi_2[[1]] = MCMC_bar_pi_2[[1]][swap,]
  
  rm(swap); rm(sel)
  
  p_vals = t(vapply(seq_len(n_genes_keep), function(gene_id){
    compute_pval_US( A = cbind(MCMC_bar_pi_1[[1]][,gene_id], MCMC_bar_pi_2[[1]][,gene_id]),
                     B = cbind(MCMC_bar_pi_1[[2]][,gene_id], MCMC_bar_pi_2[[2]][,gene_id]))
  }, FUN.VALUE = numeric(10)))
  
  rm(MCMC_bar_pi_1); rm(MCMC_bar_pi_2)
  
  p_adj = p.adjust(p_vals[,1], method = "BH")
  
  # ADD gene id for completeness!
  RES = data.frame(Transcript_id = tr_ids_keep,
                   p_val = p_vals[,1],
                   p_adj = p_adj,
                   p_vals[,-1])
  
  colnames(RES)[-c(seq_len(3))] = c("Prob-gr_B-UP",
                                    "pi_S-gr_A",
                                    "pi_U-gr_A",
                                    "pi_S-gr_B",
                                    "pi_U-gr_B",
                                    "sd_S-gr_A",
                                    "sd_U-gr_A",
                                    "sd_S-gr_B",
                                    "sd_U-gr_B")
  
  # Replace group names, with those provided in design:
  names = colnames(RES)
  sel_A = grep("gr_A", names )
  sel_B = grep("gr_B", names )
  
  colnames(RES)[sel_A] = gsub("gr_A", levels_groups[1], names[sel_A] )
  colnames(RES)[sel_B] = gsub("gr_B", levels_groups[2], names[sel_B] )
  
  # order results by significance (raw p-value)
  ord = order(RES[,4]) # 4-th column = Prob-gr_B-UP
  if(! ( any(is.na(ord)) | any(is.null(ord)) | any(is.nan(ord)) ) ){
    RES = RES[ ord, ]
  }
  
  if(traceplot){
    res = list( Differential_results = RES[, seq_len(12)],
                Convergence_results = DF_convergence,
                MCMC_U = MCMC_U)
  }else{
    res = list( Differential_results = RES[, seq_len(12)],
                Convergence_results = DF_convergence)
  }
  
  res
}
