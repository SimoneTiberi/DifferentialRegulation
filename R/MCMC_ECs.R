MCMC_ECs = function(PB_data_prepared,
                    min_counts_per_gene_per_group,
                    N_MCMC,
                    burn_in,
                    min_counts_ECs,
                    n_samples,
                    n_samples_per_group,
                    numeric_groups,
                    genes,
                    cluster,
                    cluster_ids_kept,
                    sample_ids_per_group,
                    n_groups,
                    gene_ids_sce,
                    n_genes,
                    list_EC_gene_id_original,
                    list_EC_USA_id_original,
                    cores_equal_clusters,
                    undersampling_int,
                    n_cores,
                    levels_groups,
                    traceplot){
  c_prop = 0.2 # proportionality constant of the ARW
  
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
                           
                           N_MCMC_one_cl = N_MCMC
                           burn_in_one_cl = burn_in
                           
                           counts = PB_data_prepared[[cl]][[1]]
                           S = PB_data_prepared[[cl]][[2]]
                           U = PB_data_prepared[[cl]][[3]]
                           A = PB_data_prepared[[cl]][[4]]
                           
                           if(cores_equal_clusters){
                             rm(PB_data_prepared)
                           }else{
                             PB_data_prepared[[cl]] = 1
                           }
                           #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                           # separate uniquely mapping ECs:
                           #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                           list_EC_gene_id = list_EC_gene_id_original
                           list_EC_USA_id = list_EC_USA_id_original
                           if(cores_equal_clusters){
                             rm(list_EC_gene_id_original)
                             rm(list_EC_USA_id_original)
                           }
                           
                           list_X_unique = list()
                           
                           # separate uniquely mapping ECs:
                           for(sample in seq_len(n_samples)){
                             # fill X with 0's (later replace with uniquely mapping reads):
                             X = matrix( 0, nrow = n_genes, ncol = 3)
                             
                             # unique bool vector indicating uniquely mapping ECs.
                             n_EC = length(list_EC_gene_id[[sample]])
                             unique = rep(FALSE, n_EC)
                             
                             SU_id = sel = c()
                             
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
                           
                           rm(matches)
                           
                           if(cores_equal_clusters){
                             rm(genes)
                           }
                           
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
                             # pi_S and pi_U prior:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # keep genes that we actually analyze  (sel_genes)
                             keep_sce = gene_ids_sce %in% sel_genes
                             
                             S = S[keep_sce,]
                             U = U[keep_sce,]
                             A = A[keep_sce,]
                             
                             gene_id_SUA = gene_ids_sce[keep_sce]
                             
                             rm(keep_sce);
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # DRIMSeq prior:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             if(cores_equal_clusters){
                               rm(gene_ids_sce)
                             }
                             
                             # remove genes with 0 in at least S, U or A
                             # check that dispersion is 0 for those cases...
                             # ...maybe enough that 2 are non-zero to obtain dispersion estimates?
                             sel_non_zeros = (rowSums(S) > 0) & (rowSums(U) > 0) & (rowSums(A) > 0)
                             S = S[sel_non_zeros,]
                             U = U[sel_non_zeros,]
                             A = A[sel_non_zeros,]
                             
                             S_U_A = rbind(S, U, A)
                             
                             gene_2_tr = data.frame(gene_id = rep(gene_id_SUA[sel_non_zeros], 3),
                                                    transcript_id = c(paste(gene_id_SUA[sel_non_zeros], "S"), 
                                                                      paste(gene_id_SUA[sel_non_zeros], "U"), 
                                                                      paste(gene_id_SUA[sel_non_zeros], "A") ))
                             rownames(S_U_A) = gene_2_tr$transcript_id
                             
                             log_precision = prior_precision(gene_to_transcript = gene_2_tr,
                                                             transcript_counts = S_U_A,
                                                             max_n_genes_used = 10^3,
                                                             n_cores = 1)[[2]]
                             
                             rm(gene_2_tr); rm(S_U_A)
                             
                             pi_S = rowSums(S)/( rowSums(S) + rowSums(U) + rowSums(A) )
                             pi_U = rowSums(U)/( rowSums(S) + rowSums(U) + rowSums(A) )
                             matches = match(gene_id_SUA[sel_non_zeros], names(log_precision))
                             
                             rm(S); rm(U); rm(A); rm(gene_id_SUA); rm(sel_non_zeros);
                             
                             log_delta_S = log( pi_S * exp(log_precision[matches]))
                             log_delta_U = log( pi_U * exp(log_precision[matches]))
                             
                             rm(pi_S); rm(pi_U); rm(matches)
                             
                             mean_prior = c(mean(log_precision, na.rm = TRUE),
                                            mean(log_delta_S, na.rm = TRUE),
                                            mean(log_delta_U, na.rm = TRUE))
                             sd_prior = c(sd(log_precision, na.rm = TRUE), 
                                          sd(log_delta_S, na.rm = TRUE),
                                          sd(log_delta_U, na.rm = TRUE))
                             
                             rm(log_precision); rm(log_delta_S); rm(log_delta_U); 
                             
                             # if NA or NULL or Inf, use vaguely informative values:
                             cond = is.na(mean_prior) | is.infinite(mean_prior) | is.null(mean_prior) | is.na(sd_prior) | is.infinite(sd_prior) | is.null(sd_prior)
                             mean_prior = ifelse(cond, c(log(3),  0,  0), mean_prior )
                             sd_prior = ifelse(cond, rep(5,3), sd_prior )
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # initialize objects:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # TODO: if N_MCMC_one_cl large -> use thinnings
                             MCMC_bar_pi_1 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                             MCMC_bar_pi_2 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                             MCMC_bar_pi_3 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                             
                             delta_SU = lapply(seq_len(n_groups), function(g){
                               ids = sample_ids_per_group[[g]] + 1
                               n = length(ids)
                               x = PI_SU[[ids[1]]][keep_genes_id,]
                               if(n > 1){
                                 for(i in seq.int(2, n, by = 1)){
                                   x = x + PI_SU[[ids[i]]][keep_genes_id,]
                                 }
                               }
                               x = x/n * exp(mean_prior[1])
                             })
                             
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # MCMC:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             message("Starting the MCMC")
                             
                             PI_gene_times_SU = sapply(seq_len(n_samples), function(id){
                               rep(PI_gene[,id], 3) * c(PI_SU[[id]][,1], PI_SU[[id]][,2],  PI_SU[[id]][,3]  )
                             })
                             #rm(PI_gene)
                             
                             list_X_unique = sapply(list_X_unique, c)
                             
                             # transform list_EC_gene_id, from 0 -> (n_tr-1) to 0 -> (2*n_tr-1)
                             # if S: keep number in EC
                             # if U: add n_genes_non_zero to number in EC
                             list_EC_gene_id = lapply(seq_len(n_samples), function(sample){
                               lapply(seq_along(list_EC_gene_id[[sample]]), function(ec){
                                 as.integer(list_EC_gene_id[[sample]][[ec]] + list_EC_USA_id[[sample]][[ec]] * n_genes_non_zero)
                               })
                             }) 
                             # list_EC_USA_id[[sample]][[ec]] * n_genes_non_zero:
                             # 0 if 0 (S)
                             # n_genes_non_zero if 1 (U)
                             # 2 * n_genes_non_zero if 2 (A)
                             
                             # gene ids as in "gene_id" file.
                             rm(list_EC_USA_id)
                             
                             sample_EC = rep(FALSE, N_MCMC_one_cl + 1)
                             sample_EC[seq.int(1, N_MCMC_one_cl, undersampling_int)] = TRUE
                             
                             sample_SU_TF = ifelse(seq_len(n_genes_non_zero) %in% keep_genes_id, 1, 0)
                             
                             X_list = lapply(seq_len(n_samples), matrix, data = 1, nrow = 2, ncol= 2)
                             
                             res = .Call(`_DifferentialRegulation_Rcpp_MCMC`,
                                         n_samples, # N samples
                                         n_genes_non_zero, # N genes_non_zero
                                         n_groups, # N groups
                                         n_genes_keep, # number of genes_non_zero to be analyzed.
                                         keep_genes_id-1, # vector indicating genes_non_zero to be analyzed (SU differential testing)
                                         numeric_groups - 1, # -1 ! # group id for every sample (must start from 0)
                                         sample_ids_per_group, # each list = vector with ids of samples 
                                         n_samples_per_group,
                                         N_MCMC_one_cl, # MCMC iter
                                         burn_in_one_cl, # burn-in
                                         PI_gene_times_SU, # prob of each gene (for every sample)
                                         PI_SU, # prob of each gene (for every sample)
                                         list_X_unique, # SU uniquely mapping counts
                                         list_EC_gene_id, # SU uniquely mapping counts
                                         counts, # EC counts (integers)
                                         MCMC_bar_pi_1,
                                         MCMC_bar_pi_2,
                                         MCMC_bar_pi_3,
                                         delta_SU,
                                         mean_prior,
                                         sd_prior,
                                         sample_EC,
                                         X_list,
                                         sample_SU_TF,
                                         c_prop) # 2 = sd_prior_non_informative in case prior_TF = FALSE
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # check convergence:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             convergence = my_heidel_diag(res, R = N_MCMC_one_cl, by. = 100, pvalue = 0.05)
                             rm(res)
                             
                             # set convergence (over-written below if it did not converge)
                             first = "converged"; second = "NOT run (1st chain converged)"
                             
                             if(convergence[1] == 0){ # if not converged, reset starting values and run a second chain (twice as long as the initial one):
                               rm(convergence)
                               
                               message("Our MCMC did not converge (according to Heidelberger and Welch's convergence diagnostic):
                                       we will now double 'N_MCMC_one_cl' and 'burn_in_one_cl', and run it a second time.")
                               
                               N_MCMC_one_cl = 2 * N_MCMC_one_cl
                               burn_in_one_cl = 2 * burn_in_one_cl
                               
                               # re-initialize objects:
                               MCMC_bar_pi_1 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                               MCMC_bar_pi_2 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                               MCMC_bar_pi_3 = lapply(seq_len(n_groups), matrix, data = 1, nrow = 2, ncol= 2)
                               
                               PI_SU = lapply(seq_len(n_samples), function(x){
                                 X = cbind(tot_spliced[,x], tot_unspliced[,x], tot_ambiguous[,x])
                                 X/rowSums(X)
                               })
                               
                               PI_gene_times_SU = sapply(seq_len(n_samples), function(id){
                                 rep(PI_gene[,id], 3) * c(PI_SU[[id]][,1], PI_SU[[id]][,2],  PI_SU[[id]][,3]  )
                               })
                               
                               rm(tot_spliced); rm(tot_unspliced); rm(tot_ambiguous); rm(PI_gene)
                               
                               delta_SU = lapply(seq_len(n_groups), function(g){
                                 ids = sample_ids_per_group[[g]] + 1
                                 n = length(ids)
                                 x = PI_SU[[ids[1]]][keep_genes_id,]
                                 if(n > 1){
                                   for(i in seq.int(2, n, by = 1)){
                                     x = x + PI_SU[[ids[i]]][keep_genes_id,]
                                   }
                                 }
                                 x = x/n * exp(mean_prior[1])
                               })
                               
                               sample_EC = rep(FALSE, N_MCMC_one_cl + 1)
                               sample_EC[seq.int(1, N_MCMC_one_cl, undersampling_int)] = TRUE
                               
                               res = .Call(`_DifferentialRegulation_Rcpp_MCMC`,
                                           n_samples, # N samples
                                           n_genes_non_zero, # N genes_non_zero
                                           n_groups, # N groups
                                           n_genes_keep, # number of genes_non_zero to be analyzed.
                                           keep_genes_id-1, # vector indicating genes_non_zero to be analyzed (SU differential testing)
                                           numeric_groups - 1, # -1 ! # group id for every sample (must start from 0)
                                           sample_ids_per_group, # each list = vector with ids of samples 
                                           n_samples_per_group,
                                           N_MCMC_one_cl, # MCMC iter
                                           burn_in_one_cl, # burn-in
                                           PI_gene_times_SU, # prob of each gene (for every sample)
                                           PI_SU, # prob of each gene (for every sample)
                                           list_X_unique, # SU uniquely mapping counts
                                           list_EC_gene_id, # SU uniquely mapping counts
                                           counts, # EC counts (integers)
                                           MCMC_bar_pi_1,
                                           MCMC_bar_pi_2,
                                           MCMC_bar_pi_3,
                                           delta_SU,
                                           mean_prior,
                                           sd_prior,
                                           sample_EC,
                                           X_list,
                                           sample_SU_TF,
                                           c_prop) # 2 = sd_prior_non_informative in case prior_TF = FALSE
                               
                               convergence = my_heidel_diag(res, R = N_MCMC_one_cl, by. = 100, pvalue = 0.05)
                               rm(res)
                               
                               if(convergence[1] == 0){ # if not converged for a 2nd time: return convergence error.
                                 first = "NOT converged"; second = "NOT converged"
                                 
                                 # create convergence DF:
                                 DF_convergence = data.frame(Cluster_id = cluster_ids_kept[cl],
                                                             burn_in = NA,
                                                             N_MCMC = N_MCMC_one_cl,
                                                             first_chain = first,
                                                             second_chain = second)
                                 
                                 message("Our algorithm did not converged, try to increase N_MCMC_one_cl.")
                                 
                                 return(list(NULL,
                                             DF_convergence))
                               }
                               
                               # if 2nd chain converged:
                               first = "NOT converged"; second = "converged"
                             }
                             
                             rm(list_EC_gene_id);
                             
                             message("MCMC completed and successfully converged.")
                             
                             # the code below, is only run if either chain has converged:
                             # increase the burn-in IF detected by "my_heidel_diag" (max burn-in = half of the chain length):
                             burn_in_one_cl = max(convergence[2]-1, burn_in_one_cl)
                             
                             # create convergence DF:
                             DF_convergence = data.frame(Cluster_id = cluster_ids_kept[cl],
                                                         burn_in = burn_in_one_cl,
                                                         N_MCMC = N_MCMC_one_cl,
                                                         first_chain = first,
                                                         second_chain = second)
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             # compute p-value:
                             #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                             if(traceplot){
                               # store results of MCMC, before swapping:
                               MCMC_U = list()
                               MCMC_U[[1]] = MCMC_bar_pi_2[[1]] + 0.5 * MCMC_bar_pi_3[[1]] 
                               MCMC_U[[2]] = MCMC_bar_pi_2[[2]] + 0.5 * MCMC_bar_pi_3[[2]] 
                               names(MCMC_U) = levels_groups[1:2]
                               MCMC_U$Gene_id = sel_genes
                               #Cluster_id = rep(cluster_ids_kept[cl], length(sel_genes))
                             }else{
                               MCMC_U = NULL
                             }
                             
                             # discard burn-in
                             sel = seq.int(from = burn_in_one_cl +1, to = N_MCMC_one_cl, by = 1)
                             MCMC_bar_pi_1[[1]] = MCMC_bar_pi_1[[1]][sel,]
                             MCMC_bar_pi_2[[1]] = MCMC_bar_pi_2[[1]][sel,]
                             MCMC_bar_pi_3[[1]] = MCMC_bar_pi_3[[1]][sel,]
                             MCMC_bar_pi_1[[2]] = MCMC_bar_pi_1[[2]][sel,]
                             MCMC_bar_pi_2[[2]] = MCMC_bar_pi_2[[2]][sel,]
                             MCMC_bar_pi_3[[2]] = MCMC_bar_pi_3[[2]][sel,]
                             
                             # swap A positions to decrease correlations:
                             R = length(sel)
                             swap = sample.int(R, R) # n indicates the nr of elements of the chain (exluded burn-in)
                             MCMC_bar_pi_1[[1]] = MCMC_bar_pi_1[[1]][swap,]
                             MCMC_bar_pi_2[[1]] = MCMC_bar_pi_2[[1]][swap,]
                             MCMC_bar_pi_3[[1]] = MCMC_bar_pi_3[[1]][swap,]
                             
                             rm(swap); rm(sel)
                             
                             p_vals = t(vapply(seq_len(n_genes_keep), function(gene_id){
                               compute_pval( A = cbind(MCMC_bar_pi_1[[1]][,gene_id], MCMC_bar_pi_2[[1]][,gene_id], MCMC_bar_pi_3[[1]][,gene_id]),
                                             B = cbind(MCMC_bar_pi_1[[2]][,gene_id], MCMC_bar_pi_2[[2]][,gene_id], MCMC_bar_pi_3[[2]][,gene_id]))
                             }, FUN.VALUE = numeric(22)))
                             
                             rm(MCMC_bar_pi_1); rm(MCMC_bar_pi_2); rm(MCMC_bar_pi_3)
                             
                             p_adj = p.adjust(p_vals[,1], method = "BH")
                             
                             RES = data.frame(Gene_id = sel_genes,
                                              Cluster_id = rep(cluster_ids_kept[cl], length(sel_genes)),
                                              p_val = p_vals[,1],
                                              p_adj.glb = NA,
                                              p_adj.loc = p_adj,
                                              p_vals[,-1])
                             
                             colnames(RES)[-c(seq_len(5))] = c("Prob-gr_B-UP",
                                                               "pi_S-gr_A",
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
                                       DF_convergence,
                                       MCMC_U))
                         }
  
  # merge results from multiple clusters, only if available
  if(length(order) > 1){
    RES = lapply(p_values_ALL, function(X){ # 1st element contains DR test
      X[[1]]
    })
    convergence_results = lapply(p_values_ALL, function(X){ # 2nd element contains convergence results
      X[[2]]
    })
    MCMC_U = lapply(p_values_ALL, function(X){ # 3rd element contains traceplots
      X[[3]]
    })
    
    rm(p_values_ALL)
    
    RES = do.call(rbind, RES)
    convergence_results = do.call(rbind, convergence_results)
  }else{
    RES = p_values_ALL[[1]][[1]]
    convergence_results = p_values_ALL[[1]][[2]]
    MCMC_U = p_values_ALL[[1]][[3]]
    
    rm(p_values_ALL)
  }
  # names of MCMC_U reflect the cell cluster:
  names(MCMC_U) = cluster_ids_kept[order]
  
  # CHECK if ALL NULL (if all clusters return null):
  if(!is.null(RES)){
    # adjust p-values GLOBALLY:
    RES$p_adj.glb = p.adjust(RES$p_val, method = "BH")
  }
  
  # create a final table of results, like in distinct.
  list(RES, convergence_results, MCMC_U)
}