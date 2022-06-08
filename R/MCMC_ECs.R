MCMC_ECs = function(PB_data_prepared,
                    min_counts_per_gene_per_group,
                    N_MCMC,
                    burn_in,
                    min_counts_ECs,
                    n_samples,
                    n_samples_per_group,
                    numeric_groups,
                    genes,
                    cl,
                    cluster_ids_kept,
                    sample_ids_per_group,
                    n_groups,
                    gene_ids_sce,
                    n_genes,
                    list_EC_gene_id,
                    list_EC_USA_id,
                    cores_equal_clusters
){
  # compute overall counts per cluster -> use this to rank highly abundant clusters first (likely more computationally intensive).
  overall_counts = vapply(PB_data_prepared, function(x){
    sum( unlist(x[[1]]) )
    #sum(sapply(x[[1]], sum)) # needed when unsing sparce vectors (unlist fails)
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
                           
                           counts = PB_data_prepared[[cl]][[1]]
                           S = PB_data_prepared[[cl]][[2]]
                           U = PB_data_prepared[[cl]][[3]]
                           A = PB_data_prepared[[cl]][[4]]
                           
                           if(cores_equal_clusters){
                             rm(PB_data_prepared)
                           }else{
                             PB_data_prepared[[cl]] = NULL
                           }
                           #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                           # separate uniquely mapping ECs:
                           #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                           list_X_unique = list()
                           
                           # separate uniquely mapping ECs:
                           for(sample in seq_len(n_samples)){
                             # fill X with 0's (later replace with uniquely mapping reads):
                             X = matrix( 0, nrow = n_genes, ncol = 3)
                             
                             # unique bool vector indicating uniquely mapping ECs.
                             n_EC = length(list_EC_gene_id[[sample]])
                             unique = rep(FALSE, n_EC)
                             
                             for(ec in seq_len(n_EC)){
                               gene_id = list_EC_gene_id[[sample]][[ec]] + 1
                               if(length(gene_id) == 1){
                                 unique[ec] = TRUE
                                 SU_id = list_EC_USA_id[[sample]][[ec]] + 1
                                 sel = cbind(gene_id, SU_id)
                                 
                                 X[ sel ] = X[ sel ] + counts[[sample]][ec]
                               }
                             }
                             
                             list_X_unique[[sample]] = X
                             
                             # remove "unique" ECs from latent variables (already assigned in "list_X_unique")
                             list_EC_gene_id[[sample]] = list_EC_gene_id[[sample]][!unique]
                             list_EC_USA_id[[sample]] = list_EC_USA_id[[sample]][!unique]
                             counts[[sample]] = counts[[sample]][!unique]
                           }
                           
                           rm(X); rm(unique); rm(sel); rm(SU_id); rm(n_EC)
                           
                           #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                           # FILTER 1) remove ECs with 0 counts (within a specific cell type):
                           #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                           # 1) remove ECs with 0 counts (within a specific cell type):
                           sel_non_zero_EC = lapply(counts, function(X){
                             X > min_counts_ECs
                           })
                           
                           # lapply(list_EC_gene_id, length); lapply(list_EC_USA_id, length); lapply(counts, length)
                           
                           # trim EC gene list, each element = 1 EC, OK:
                           list_EC_gene_id = lapply(seq_len(n_samples), function(X){
                             list_EC_gene_id[[X]][ sel_non_zero_EC[[X]] ]
                           })
                           
                           # trim EC SU list, each element = 1 EC, OK:
                           list_EC_USA_id = lapply(seq_len(n_samples), function(X){
                             list_EC_USA_id[[X]][ sel_non_zero_EC[[X]] ]
                           })
                           
                           # trim counts, each element = 1 EC, OK:
                           counts = lapply(seq_len(n_samples), function(X){
                             counts[[X]][ sel_non_zero_EC[[X]] ]
                           })
                           
                           #lapply(list_EC_gene_id, length); lapply(list_EC_USA_id, length); lapply(counts, length)
                           
                           rm(sel_non_zero_EC)
                           
                           #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                           # FILTER 2) remove genes absent across ALL ECs of ALL genes...check sce of those genes -> it'd be ~ 0!
                           #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                           # Select genes in ECs with non-zero counts:
                           genes_non_zero = genes[ unique(unlist(list_EC_gene_id)) + 1 ]
                           # PLUS ADD GENES with UNIQUELY MAPPING READS!!!
                           genes_non_zero = c(genes_non_zero,  genes[ rowSums(do.call(cbind, list_X_unique)) > 0 ] )
                           genes_non_zero = unique(genes_non_zero)
                           genes_non_zero = sort(genes_non_zero)
                           
                           # for testing purposed: line 356
                           # genes_non_zero = sel_genes
                           
                           # length(genes_non_zero)/n_genes
                           
                           # list_X_unique -> remove rows for non-selected genes
                           # list_EC_gene_id -> update gene ids
                           # n_genes -> n_genes_non_zero
                           
                           # genes
                           matches = match(genes, genes_non_zero)
                           
                           list_EC_gene_id = lapply(list_EC_gene_id, function(X){
                             lapply(X, function(x){
                               matches[x + 1] - 1
                             })
                           })
                           
                           matches = match(genes_non_zero, genes)
                           list_X_unique = lapply(list_X_unique, function(x){
                             # need to do: match(genes, genes_non_zero)
                             # genes_non_zero don't necessarily follow the same order of genes!!
                             x[matches,] #
                           })
                           
                           rm(matches); rm(genes)
                           
                           n_genes_non_zero = length(genes_non_zero)
                           
                           #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                           # set starting values, defined from "sce" values
                           #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                           TOT_counts = S + U + A
                           # gene_ids_sce is the row name of TOT_counts
                           
                           # use pi_gene, estimated from sample-specific counts (+ 1):
                           matches = match(genes_non_zero, gene_ids_sce)
                           
                           PI_gene = 1 + TOT_counts[matches, ]
                           rownames(PI_gene) = genes_non_zero
                           PI_gene[is.na(PI_gene)] = 1 # set NAs to 1
                           PI_gene[ PI_gene < 1 ] = 1 # set counts < 1 to 1
                           PI_gene = apply(PI_gene, 2, function(x) x/sum(x))
                           
                           # assign pi_SU as starting MCMC value:
                           tot_spliced = 1/3 + S
                           tot_spliced = tot_spliced[matches, ]
                           tot_spliced[is.na(tot_spliced)] = 1/3
                           
                           tot_unspliced = 1/3 + U
                           tot_unspliced = tot_unspliced[matches, ]
                           tot_unspliced[is.na(tot_unspliced)] = 1/3
                           
                           tot_ambiguous = 1/3 + A
                           tot_ambiguous = tot_ambiguous[matches, ]
                           tot_ambiguous[is.na(tot_ambiguous)] = 1/3
                           
                           PI_SU = lapply(seq_len(n_samples), function(x){
                             X = cbind(tot_spliced[,x], tot_unspliced[,x], tot_ambiguous[,x])
                             X/rowSums(X)
                           })
                           
                           rm(matches)
                           
                           # set filter to analyze genes_non_zero: at least xxx counts across all cells.
                           counts_per_group = vapply(sample_ids_per_group, function(id){
                             rowSums(TOT_counts[, id + 1])
                           }, FUN.VALUE = numeric( nrow(TOT_counts) ) )
                           
                           rm(TOT_counts)
                           
                           sel_genes = rowSums(counts_per_group >= min_counts_per_gene_per_group) == n_groups
                           sel_genes = gene_ids_sce[sel_genes]
                           
                           rm(counts_per_group)
                           
                           # to guarantee that the order is preserved:
                           sel_genes = genes_non_zero[ genes_non_zero %in% sel_genes ]
                           
                           keep_genes_id = which(genes_non_zero %in% sel_genes)
                           n_genes_keep = length(keep_genes_id)
                           n_genes_keep
                           
                           rm(genes_non_zero)
                           #n_genes_keep; n_genes_non_zero
                           
                           if(n_genes_keep > 0){ # if at least 1 gene is selected:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # DRIMSeq prior:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # Infer pseudo-bulk counts:
                             sel_genes_random = sample(sel_genes, min(10^2, length(sel_genes)), replace = FALSE)
                             
                             keep_sce = gene_ids_sce %in% sel_genes_random
                             
                             S = S[keep_sce,]
                             U = U[keep_sce,]
                             A = A[keep_sce,]
                             
                             #summary(rowSums(S + U + A))
                             # should be > min selected above!!!
                             
                             # compute prior for the dispersion via DRIMSeq, based on the selected genes_non_zero only!
                             
                             gene_id_SUA = gene_ids_sce[keep_sce]
                             
                             rm(keep_sce); rm(sel_genes_random); rm(gene_ids_sce)
                             
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
                             # only store mean and var of precision.
                             
                             rm(gene_2_tr); rm(S_U_A)
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # initialize objects:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # TODO: if N_MCMC large -> use thinnings
                             MCMC_bar_pi_1 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                             MCMC_bar_pi_2 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                             MCMC_bar_pi_3 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                             
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
                             
                             res = .Call(`_DifferentialRegulation_Rcpp_MCMC`,
                                         n_samples, # N samples
                                         n_genes_non_zero, # N genes_non_zero
                                         n_groups, # N groups
                                         n_genes_keep, # number of genes_non_zero to be analyzed.
                                         keep_genes_id-1, # vector indicating genes_non_zero to be analyzed (SU differential testing)
                                         numeric_groups - 1, # -1 ! # group id for every sample (must start from 0)
                                         sample_ids_per_group, # each list = vector with ids of samples 
                                         n_samples_per_group,
                                         N_MCMC, # MCMC iter
                                         burn_in, # burn-in
                                         PI_gene, # prob of each gene (for every sample)
                                         PI_SU, # prob of each gene (for every sample)
                                         list_X_unique, # SU uniquely mapping counts
                                         list_EC_gene_id, # SU uniquely mapping counts
                                         list_EC_USA_id, # TRUE -> S; FALSE -> U
                                         counts, # EC counts (integers)
                                         MCMC_bar_pi_1,
                                         MCMC_bar_pi_2,
                                         MCMC_bar_pi_3,
                                         chol,
                                         delta_SU,
                                         TRUE, # I ALWAYS USE THE PRIOR
                                         precision[1],
                                         precision[2], 
                                         2) # 2 = sd_prior_non_informative in case prior_TF = FALSE
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # check convergence:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             convergence = my_heidel_diag(res, R = N_MCMC, by. = 100, pvalue = 0.01)
                             rm(res)
                             
                             # set convergence (over-written below if it did not converge)
                             first = "converged"; second = NULL   
                             
                             if(convergence[1] == 0){ # if not converged, reset starting values and run a second chain (twice as long as the initial one):
                               rm(convergence)
                               
                               message("Our MCMC did not converge (according to Heidelberger and Welch's convergence diagnostic):
                                       we will now double 'N_MCMC' and 'burn_in', and run it a second time.")
                               
                               N_MCMC = 2 * N_MCMC
                               burn_in = 2 * burn_in
                               
                               # re-initialize objects:
                               MCMC_bar_pi_1 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                               MCMC_bar_pi_2 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                               MCMC_bar_pi_3 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                               
                               PI_SU = lapply(seq_len(n_samples), function(x){
                                 X = cbind(tot_spliced[,x], tot_unspliced[,x], tot_ambiguous[,x])
                                 X/rowSums(X)
                               })
                               
                               rm(tot_spliced); rm(tot_unspliced); rm(tot_ambiguous)
                               
                               delta_SU = lapply(seq_len(n_groups), function(g){
                                 ids = sample_ids_per_group[[g]] + 1
                                 n = length(ids)
                                 x = PI_SU[[ids[1]]]
                                 if(n > 1){
                                   for(i in seq.int(2, n, by = 1)){
                                     x = x + PI_SU[[ids[i]]]
                                   }
                                 }
                                 x = x/n * exp(precision$prior[1])
                               })
                               chol = lapply(seq_len(n_groups), function(i){
                                 lapply(seq_len(n_genes_keep), matrix, data = 1, nrow = 3, ncol= 3)
                               })
                               
                               res = .Call(`_DifferentialRegulation_Rcpp_MCMC`,
                                           n_samples, # N samples
                                           n_genes_non_zero, # N genes_non_zero
                                           n_groups, # N groups
                                           n_genes_keep, # number of genes_non_zero to be analyzed.
                                           keep_genes_id-1, # vector indicating genes_non_zero to be analyzed (SU differential testing)
                                           numeric_groups - 1, # -1 ! # group id for every sample (must start from 0)
                                           sample_ids_per_group, # each list = vector with ids of samples 
                                           n_samples_per_group,
                                           N_MCMC, # MCMC iter
                                           burn_in, # burn-in
                                           PI_gene, # prob of each gene (for every sample)
                                           PI_SU, # prob of each gene (for every sample)
                                           list_X_unique, # SU uniquely mapping counts
                                           list_EC_gene_id, # SU uniquely mapping counts
                                           list_EC_USA_id, # TRUE -> S; FALSE -> U
                                           counts, # EC counts (integers)
                                           MCMC_bar_pi_1,
                                           MCMC_bar_pi_2,
                                           MCMC_bar_pi_3,
                                           chol,
                                           delta_SU,
                                           TRUE, # I ALWAYS USE THE PRIOR
                                           precision$prior[1],
                                           precision$prior[2], 
                                           2) # 2 = sd_prior_non_informative in case prior_TF = FALSE
                               
                               convergence = my_heidel_diag(res, R = N_MCMC, by. = 100, pvalue = 0.01)
                               rm(res)
                               
                               if(convergence[1] == 0){ # if not converged for a 2nd time: return convergence error.
                                 first = "NOT converged"; second = "NOT converged"
                                 
                                 # create convergence DF:
                                 DF_convergence = data.frame(Cluster_id = cluster_ids_kept[cl],
                                                             burn_in = NA,
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
                             
                             message("MCMC completed and successfully converged.")
                             
                             # the code below, is only run if either chain has converged:
                             # increase the burn-in IF detected by "my_heidel_diag" (max burn-in = half of the chain length):
                             burn_in = max(convergence[2]-1, burn_in)
                             
                             # create convergence DF:
                             DF_convergence = data.frame(Cluster_id = cluster_ids_kept[cl],
                                                         burn_in = burn_in,
                                                         N_MCMC = N_MCMC,
                                                         first_chain = first,
                                                         second_chain = second)
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # compute p-value:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             sel = seq.int(from = burn_in +1, to = N_MCMC, by = 1)
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
    RES = p_values_ALL[[1]]
    convergence_results = p_values_ALL[[2]]
  }
  
  # CHECK if ALL NULL (if all clusters return null):
  if(!is.null(RES)){
    # adjust p-values GLOBALLY:
    RES$p_adj.glb = p.adjust(RES$p_val, method = "BH")
  }
  
  # create a final table of results, like in distinct.
  RES
}