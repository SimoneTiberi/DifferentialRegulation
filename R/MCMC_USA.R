MCMC_USA = function(PB_data_prepared,
                    min_counts_per_gene_per_group,
                    N_MCMC,
                    burn_in,
                    n_samples,
                    n_samples_per_group,
                    numeric_groups,
                    cluster,
                    cluster_ids_kept,
                    sample_ids_per_group,
                    n_groups,
                    gene_ids_sce,
                    cores_equal_clusters
){
  # compute overall counts per cluster -> use this to rank highly abundant clusters first (likely more computationally intensive).
  overall_counts = vapply(PB_data_prepared, function(x){
    sum( unlist(x[ seq.int(2, 4, by = 1) ]) )
  }, FUN.VALUE = numeric(1))
  order = order(overall_counts, decreasing = TRUE)
  
  rm(overall_counts)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # run in Parallel:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  p_values_ALL = foreach(cl = order,
                         # maybe only via R package
                         .packages=c("DifferentialRegulation"),
                         .errorhandling = "stop") %dorng%{
                           
                           N_MCMC_one_cl = N_MCMC
                           burn_in_one_cl = burn_in
                           
                           S = PB_data_prepared[[cl]][[2]]
                           U = PB_data_prepared[[cl]][[3]]
                           A = PB_data_prepared[[cl]][[4]]
                           
                           if(cores_equal_clusters){
                             rm(PB_data_prepared)
                           }else{
                             PB_data_prepared[[cl]] = 1
                           }
                           
                           SUA = list()
                           for(i in seq_len(n_samples)){
                             SUA[[i]] = cbind(S[,i], U[,i], A[,i])
                           }
                           
                           n_genes = nrow(SUA[[1]])
                           
                           # set filter to analyze genes_non_zero: at least xxx counts across all cells.
                           sample_counts = vapply(SUA, rowSums, FUN.VALUE = numeric( n_genes ) )
                           counts_per_group = vapply(sample_ids_per_group, function(id){
                             rowSums(sample_counts[,id + 1])
                           }, FUN.VALUE = numeric( n_genes ) )
                           
                           rm(sample_counts)
                           
                           sel = rowSums(counts_per_group >= min_counts_per_gene_per_group) == n_groups
                           sel_genes = gene_ids_sce[sel]
                           
                           rm(counts_per_group)
                           
                           # to guarantee that the order is preserved:
                           sel_genes = gene_ids_sce[ gene_ids_sce %in% sel_genes ]
                           
                           n_genes_keep = length(sel_genes)
                           
                           if(n_genes_keep > 0){ # if at least 1 gene is selected:
                             # filter SUA object according to genes that pass filtering:
                             for(i in seq_len(n_samples)){
                               SUA[[i]] = SUA[[i]][sel,]
                             }
                             S = S[sel,]
                             U = U[sel,]
                             A = A[sel,]
                             
                             rm(sel)
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # DRIMSeq prior:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # Infer pseudo-bulk counts:
                             sel_genes_random = sample(sel_genes, min(10^2, length(sel_genes)), replace = FALSE)
                             
                             keep_sce = sel_genes %in% sel_genes_random
                             
                             rm(sel_genes_random)
                             
                             S = S[keep_sce,]
                             U = U[keep_sce,]
                             A = A[keep_sce,]
                             
                             gene_id_SUA = sel_genes[keep_sce]
                             rm(keep_sce)
                             
                             S_U_A = rbind(S, U, A)
                             
                             rm(S); rm(U); rm(A);
                             
                             gene_2_tr = data.frame(gene_id = rep(gene_id_SUA, 3),
                                                    transcript_id = c(paste(gene_id_SUA, "S"), 
                                                                      paste(gene_id_SUA, "U"), 
                                                                      paste(gene_id_SUA, "A") ))
                             rownames(S_U_A) = gene_2_tr$transcript_id
                             
                             precision = prior_precision(gene_to_transcript = gene_2_tr,
                                                         transcript_counts = S_U_A,
                                                         n_cores = 1)[[1]]
                             
                             rm(gene_2_tr); rm(S_U_A)
                             
                             # if NA or NULL or Inf, use vaguely informative values:
                             if(is.na(precision[1]) | is.infinite(precision[1]) | is.null(precision[1])){
                               precision[1] = 3
                             }
                             if(is.na(precision[2]) | is.infinite(precision[2]) | is.null(precision[2])){
                               precision[2] = 10
                             }
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # initialize objects:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             MCMC_bar_pi_1 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                             MCMC_bar_pi_2 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                             MCMC_bar_pi_3 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                             
                             # assign pi_SU as starting MCMC value:
                             PI_SU = lapply(seq_len(n_samples), function(x){
                               X = 1/3 + SUA[[i]]
                               X/rowSums(X)
                             })
                             
                             delta_SU = lapply(seq_len(n_groups), function(g){
                               ids = sample_ids_per_group[[g]] + 1
                               n = length(ids)
                               x = PI_SU[[ids[1]]]
                               if(n > 1){
                                 for(i in seq.int(2, n, by = 1)){
                                   x = x + PI_SU[[ids[i]]]
                                 }
                               }
                               x = x/n * exp(precision[1])
                             })
                             
                             chol = lapply(seq_len(n_groups), function(i){
                               lapply(seq_len(n_genes_keep), matrix, data = 1, nrow = 3, ncol= 3)
                             })
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # MCMC:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             message("Starting the MCMC")
                             
                             res = .Call(`_DifferentialRegulation_Rcpp_MCMC_sce`,
                                         n_samples, # N samples
                                         n_genes_keep, # N genes_non_zero
                                         n_groups, # N groups
                                         numeric_groups - 1, # -1 ! # groups id for every sample (must start from 0)
                                         sample_ids_per_group, # each list = vector with ids of samples 
                                         n_samples_per_group,
                                         N_MCMC_one_cl, # MCMC iter
                                         burn_in_one_cl, # burn-in
                                         PI_SU, # prob of each gene (for every sample)
                                         SUA, # SU uniquely mapping counts
                                         MCMC_bar_pi_1,
                                         MCMC_bar_pi_2,
                                         MCMC_bar_pi_3,
                                         chol,
                                         delta_SU,
                                         TRUE, # I ALWAYS USE THE PRIOR
                                         precision[1],
                                         precision[2], 
                                         2)
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # check convergence:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             convergence = my_heidel_diag(res, R = N_MCMC_one_cl, by. = 100, pvalue = 0.01)
                             rm(res)
                             
                             # set convergence (over-written below if it did not converge)
                             first = "converged"; second = "NOT run (1st chain converged)"
                             
                             if(convergence[1] == 0){ # if not converged, reset starting values and run a second chain (twice as long as the initial one):
                               message("Our MCMC did not converge (according to Heidelberger and Welch's convergence diagnostic):
                                       we will now double 'N_MCMC_one_cl' and 'burn_in_one_cl', and run it a second time.")
                               
                               N_MCMC_one_cl = 2 * N_MCMC_one_cl
                               burn_in_one_cl = 2 * burn_in_one_cl
                               
                               # re-initialize objects:
                               MCMC_bar_pi_1 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                               MCMC_bar_pi_2 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                               MCMC_bar_pi_3 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                               
                               # assign pi_SU as starting MCMC value:
                               PI_SU = lapply(seq_len(n_samples), function(x){
                                 X = 1/3 + SUA[[i]]
                                 X/rowSums(X)
                               })
                               
                               delta_SU = lapply(seq_len(n_groups), function(g){
                                 ids = sample_ids_per_group[[g]] + 1
                                 n = length(ids)
                                 x = PI_SU[[ids[1]]]
                                 if(n > 1){
                                   for(i in seq.int(2, n, by = 1)){
                                     x = x + PI_SU[[ids[i]]]
                                   }
                                 }
                                 x = x/n * exp(precision[1])
                               })
                               
                               chol = lapply(seq_len(n_groups), function(i){
                                 lapply(seq_len(n_genes_keep), matrix, data = 1, nrow = 3, ncol= 3)
                               })
                               
                               res = .Call(`_DifferentialRegulation_Rcpp_MCMC_sce`,
                                           n_samples, # N samples
                                           n_genes_keep, # N genes_non_zero
                                           n_groups, # N groups
                                           numeric_groups - 1, # -1 ! # groups id for every sample (must start from 0)
                                           sample_ids_per_group, # each list = vector with ids of samples 
                                           n_samples_per_group,
                                           N_MCMC_one_cl, # MCMC iter
                                           burn_in_one_cl, # burn-in
                                           PI_SU, # prob of each gene (for every sample)
                                           SUA, # SU uniquely mapping counts
                                           MCMC_bar_pi_1,
                                           MCMC_bar_pi_2,
                                           MCMC_bar_pi_3,
                                           chol,
                                           delta_SU,
                                           TRUE, # I ALWAYS USE THE PRIOR
                                           precision[1],
                                           precision[2], 
                                           2)
                               
                               convergence = my_heidel_diag(res, R = N_MCMC_one_cl, by. = 100, pvalue = 0.01)
                               rm(res)
                               
                               if(convergence[1] == 0){ # if not converged for a 2nd time: return convergence error.
                                 first = "NOT converged"; second = "NOT converged"
                                 
                                 # create convergence DF:
                                 DF_convergence = data.frame(Cluster_id = cluster_ids_kept[cl],
                                                             burn_in_one_cl = NA,
                                                             N_MCMC_one_cl = N_MCMC_one_cl,
                                                             first_chain = first,
                                                             second_chain = second)
                                 
                                 message("Our algorithm did not converged, try to increase N_MCMC_one_cl.")
                                 
                                 return(list(NULL,
                                             DF_convergence))
                               }
                               
                               # if 2nd chain converged:
                               first = "NOT converged"; second = "converged"
                             }
                             
                             message("MCMC completed and successfully converged.")
                             
                             # the code below, is only run if either chain has converged:
                             # increase the burn-in IF detected by "my_heidel_diag" (max burn-in = half of the chain length):
                             burn_in_one_cl = max(convergence[2]-1, burn_in_one_cl)
                             
                             # create convergence DF:
                             DF_convergence = data.frame(Cluster_id = cluster_ids_kept[cl],
                                                         burn_in_one_cl = burn_in_one_cl,
                                                         N_MCMC_one_cl = N_MCMC_one_cl,
                                                         first_chain = first,
                                                         second_chain = second)
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # compute p-value:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             sel = seq.int(from = burn_in_one_cl +1, to = N_MCMC_one_cl, by = 1)
                             p_vals = t(vapply(seq_len(n_genes_keep), function(gene_id){
                               a = MCMC_bar_pi_1[[1]][sel,gene_id]
                               b = MCMC_bar_pi_2[[1]][sel,gene_id]
                               c = MCMC_bar_pi_3[[1]][sel,gene_id]
                               tot = a+b+c
                               A = cbind(a, b, c)/tot
                               
                               a = MCMC_bar_pi_1[[2]][sel,gene_id]
                               b = MCMC_bar_pi_2[[2]][sel,gene_id]
                               c = MCMC_bar_pi_3[[2]][sel,gene_id]
                               tot = a+b+c
                               B = cbind(a, b, c)/tot
                               
                               compute_pval( A = A, B = B, K = 3, N = n_samples)
                             }, FUN.VALUE = numeric(21)))
                             # TODO: speed-up p-val computation!
                             
                             rm(MCMC_bar_pi_1); rm(MCMC_bar_pi_2); rm(MCMC_bar_pi_3)
                             
                             p_adj = p.adjust(p_vals[,1], method = "BH")
                             
                             RES = data.frame(Gene_id = sel_genes,
                                              Cluster_id = rep(cluster_ids_kept[cl], length(sel_genes)),
                                              p_val = p_vals[,1],
                                              p_adj.glb = NA,
                                              p_adj.loc = p_adj,
                                              p_vals[,-1])
                             
                             colnames(RES)[-c(seq_len(5))] = c("pi_S-gr_A",
                                                               "pi_U-gr_A",
                                                               "pi_S-gr_B",
                                                               "pi_U-gr_B",
                                                               "sd_S-gr_A",
                                                               "sd_U-gr_A",
                                                               "sd_S-gr_B",
                                                               "sd_U-gr_B",
                                                               "pi_S-gr_A",
                                                               "pi_U-gr_A",
                                                               "pi_A-gr_A",
                                                               "pi_S-gr_B",
                                                               "pi_U-gr_B",
                                                               "pi_A-gr_B",
                                                               "sd_S-gr_A",
                                                               "sd_U-gr_A",
                                                               "sd_A-gr_A",
                                                               "sd_S-gr_B",
                                                               "sd_U-gr_B",
                                                               "sd_A-gr_B")
                           }else{ # if n_genes_keep == 0
                             RES = NULL
                             DF_convergence = NULL
                           }
                           return(list(RES,
                                       DF_convergence))
                         }
  
  # merge results from multiple clusters, only if available
  if(length(order) > 1){
    p_values_ALL_test = lapply(p_values_ALL, function(X){ # 1st element contains DR test
      X[[1]]
    })
    p_values_ALL_convergence = lapply(p_values_ALL, function(X){ # 2nd element contains convergence results
      X[[2]]
    })
    
    RES = do.call(rbind, p_values_ALL_test)
    
    convergence_results = do.call(rbind, p_values_ALL_convergence)
  }else{
    RES = p_values_ALL[[1]][[1]]
    convergence_results = p_values_ALL[[1]][[2]]
  }
  
  # CHECK if ALL NULL (if all clusters return null):
  if(!is.null(RES)){
    # adjust p-values GLOBALLY:
    RES$p_adj.glb = p.adjust(RES$p_val, method = "BH")
  }
  
  list(RES, convergence_results)
}