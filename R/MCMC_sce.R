MCMC_sce = function(PB_data_prepared,
                    samples_design,
                    min_counts_per_gene_per_group,
                    N_MCMC,
                    burn_in,
                    n_samples,
                    n_samples_per_group,
                    numeric_groups,
                    genes,
                    cl,
                    cluster_ids_kept,
                    sample_ids_per_group,
                    n_groups,
                    gene_ids_sce
){
  # compute overall counts per cluster -> use this to rank highly abundant clusters first (likely more computationally intensive).
  overall_counts = sapply(PB_data_prepared, function(x){
    sum( unlist(x[[2:4]]) )
  })
  order = order(overall_counts, decreasing = TRUE)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # run in Parallel:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  p_values_ALL = foreach(cl = order,
                         # maybe only via R package
                         .packages=c("DifferentialRegulation"),
                         .errorhandling = "stop") %dorng%{
                           
                           S = PB_data_prepared[[cl]][[2]]
                           U = PB_data_prepared[[cl]][[3]]
                           A = PB_data_prepared[[cl]][[4]]
                           
                           SUA = list()
                           for(i in 1:n_samples){
                             SUA[[i]] = cbind(S[,i], U[,i], A[,i])
                           }
                           
                           # set filter to analyze genes_non_zero: at least xxx counts across all cells.
                           sample_counts = sapply(SUA, rowSums)
                           counts_per_group = sapply(sample_ids_per_group, function(id){
                             rowSums(sample_counts[,id + 1])
                           })
                           
                           sel = rowSums(counts_per_group >= min_counts_per_gene_per_group) == n_groups
                           sel_genes = gene_ids_sce[sel_genes]
                           
                           # to guarantee that the order is preserved:
                           sel_genes = genes[ genes %in% sel_genes ]
                           
                           n_genes_keep = length(sel_genes)
                           
                           message(paste("n_genes_keep:", n_genes_keep)) 
                           
                           if(n_genes_keep > 0){ # if at least 1 gene is selected:
                             # filter SUA object according to genes that pass filtering:
                             for(i in 1:n_samples){
                               SUA[[i]] = SUA[[i]][sel,]
                               S = S[sel,]
                               U = U[sel,]
                               A = A[sel,]
                             }
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # DRIMSeq prior:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # Infer pseudo-bulk counts:
                             sel_genes_random = sample(sel_genes, min(10^2, length(sel_genes)), replace = FALSE)
                             
                             keep_sce = genes %in% sel_genes_random
                             
                             S = S[keep_sce,]
                             U = U[keep_sce,]
                             A = A[keep_sce,]
                             
                             gene_id_SUA = genes[keep_sce]
                             
                             S_U_A = rbind(S, U, A)
                             gene_2_tr = data.frame(gene_id = rep(gene_id_SUA, 3),
                                                    transcript_id = c(paste(gene_id_SUA, "S"), 
                                                                      paste(gene_id_SUA, "U"), 
                                                                      paste(gene_id_SUA, "A") ))
                             rownames(S_U_A) = gene_2_tr$transcript_id
                             
                             precision = prior_precision(gene_to_transcript = gene_2_tr,
                                                         transcript_counts = S_U_A,
                                                         n_cores = 1)
                             rm(keep_sce)
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # initialize objects:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             thinning = as.integer(1)
                             
                             MCMC_bar_pi_1 = lapply(1:n_groups, matrix, data = 1, nrow = 2, ncol= 2)
                             MCMC_bar_pi_2 = lapply(1:n_groups, matrix, data = 1, nrow = 2, ncol= 2)
                             MCMC_bar_pi_3 = lapply(1:n_groups, matrix, data = 1, nrow = 2, ncol= 2)
                             
                             # assign pi_SU as starting MCMC value:
                             PI_SU = lapply(1:n_samples, function(x){
                               X = 1/3 + SUA[[i]]
                               X/rowSums(X)
                             })
                             
                             delta_SU = lapply(1:n_groups, function(g){
                               ids = sample_ids_per_group[[g]] + 1
                               n = length(ids)
                               x = PI_SU[[ids[1]]]
                               if(n > 1){
                                 for(i in 2:n){
                                   x = x + PI_SU[[ids[i]]]
                                 }
                               }
                               x = x/n * exp(precision$prior[1])
                             })
                             
                             chol = lapply(1:n_groups, function(i){
                               lapply(1:n_genes, matrix, data = 1, nrow = 3, ncol= 3)
                             })
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # MCMC:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             res = MCMC_no_ECs( n_samples, # N samples
                                                n_genes, # N genes_non_zero
                                                n_groups, # N groups
                                                numeric_groups - 1, # -1 ! # groups id for every sample (must start from 0)
                                                sample_ids_per_group, # each list = vector with ids of samples 
                                                n_samples_per_group,
                                                N_MCMC, # MCMC iter
                                                burn_in, # burn-in
                                                thinning,
                                                PI_SU, # prob of each gene (for every sample)
                                                SUA, # SU uniquely mapping counts
                                                MCMC_bar_pi_1,
                                                MCMC_bar_pi_2,
                                                MCMC_bar_pi_3,
                                                chol,
                                                delta_SU,
                                                prior_TF,
                                                precision$prior[1],
                                                precision$prior[2], 
                                                2)
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # check convergence:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             convergence = my_heidel_diag(res[[1]], R = N_MCMC, by. = 100, pvalue = 0.01)
                             
                             if(convergence[1] == 0){ # if not converged, reset starting values and run a second chain (twice as long as the initial one):
                               N_MCMC = 2 * N_MCMC
                               burn_in = 2 * burn_in
                               
                               PI_SU = lapply(1:n_samples, function(i){
                                 x = SUA[[i]] + 1
                                 matrix(x/rowSums(x), ncol = 3, byrow = FALSE)
                               })
                               delta_SU = lapply(1:n_groups, function(g){
                                 ids = sample_ids_per_group[[g]] + 1
                                 n = length(ids)
                                 x = PI_SU[[ids[1]]]
                                 if(n > 1){
                                   for(i in 2:n){
                                     x = x + PI_SU[[ids[i]]]
                                   }
                                 }
                                 x = x/n * exp(precision$prior[1])
                               })
                               
                               res = MCMC_no_ECs( n_samples, # N samples
                                                  n_genes, # N genes_non_zero
                                                  n_groups, # N groups
                                                  numeric_groups - 1, # -1 ! # groups id for every sample (must start from 0)
                                                  sample_ids_per_group, # each list = vector with ids of samples 
                                                  n_samples_per_group,
                                                  N_MCMC, # MCMC iter
                                                  burn_in, # burn-in
                                                  thinning,
                                                  PI_SU, # prob of each gene (for every sample)
                                                  SUA, # SU uniquely mapping counts
                                                  MCMC_bar_pi_1,
                                                  MCMC_bar_pi_2,
                                                  MCMC_bar_pi_3,
                                                  chol,
                                                  delta_SU,
                                                  prior_TF,
                                                  precision$prior[1],
                                                  precision$prior[2], 
                                                  2)
                               
                               convergence = my_heidel_diag(res[[1]], R = N_MCMC, by. = 100, pvalue = 0.01)
                               
                               if(convergence[1] == 0){ # if not converged for a 2nd time: return convergence error.
                                 return("Our algorithm did not converged, try to increase N_MCMC.")
                               }
                             }
                             # the code below, is only run if either chain has converged:
                             # increase the burn-in IF detected by "my_heidel_diag" (max burn-in = half of the chain length):
                             burn_in = max(convergence[2]-1, burn_in)
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # compute p-value:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             p_vals = sapply(1:n_genes_keep, function(gene_id){
                               a = res[[2]][[1]][(burn_in+1):N_MCMC,gene_id]
                               b = res[[3]][[1]][(burn_in+1):N_MCMC,gene_id]
                               c = res[[4]][[1]][(burn_in+1):N_MCMC,gene_id]
                               tot = a+b+c
                               A = cbind(a, b, c)/tot
                               
                               a = res[[2]][[2]][(burn_in+1):N_MCMC,gene_id]
                               b = res[[3]][[2]][(burn_in+1):N_MCMC,gene_id]
                               c = res[[4]][[2]][(burn_in+1):N_MCMC,gene_id]
                               tot = a+b+c
                               B = cbind(a, b, c)/tot
                               
                               compute_pval_FULL( A = A, B = B, K = 3, N = n_samples)
                             })
                             p_adj = p.adjust(p_vals, method = "BH")
                             
                             RES = data.frame(Gene_id = sel_genes,
                                              Cluster_id = rep(cluster_ids_kept[cl], length(sel_genes)),
                                              p_val = p_vals,
                                              p_adj.loc = p_adj)
                           }else{ # if n_genes_keep == 0
                             RES = NULL
                           }
                           return(RES)
                         }
  
  RES = do.call(rbind, p_values_ALL)
  
  # adjust p-values GLOBALLY:
  RES$p_adj.glb = p.adjust(RES$p_vals, method = "BH")
  
  # create a final table of results, like in distinct.
  RES
}