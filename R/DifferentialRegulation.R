#' Discover differentially regulated genes
#'
#' \code{DifferentialRegulation} identified differentially regulated genes between two conditions 
#' (e.g., healthy vs. disease or treated vs. untreated) in each cluster of cells.
#' Parameters are inferred via Markov chain Monte Carlo (MCMC) techniques and a differential testing is performed 
#' via a multivariate Wald test on the posterior densities of the group-level USA (Unspliced, Spliced and Ambiguous) counts relative abundance.
#' 
#' @param PB_counts a \code{list}, computed via \code{\link{compute_PB_counts}}
#' @param EC a logical, indicating whether to use equivalence classes (if TRUE, default)
#' or USA estimated counts (if FALSE).
#' @param n_cores the number of cores to parallelize the tasks on.
#' Since parallelization is at the cluster level (each cluster is parallelized on a thread), 
#' we suggest setting n_cores to the number of clusters (e.g., cell-types), as set by default if 'n_cores' is not specified.
#' @param N_MCMC the number of iterations for the MCMC algorithm (including burn-in). Min 2*10^3.
#' If our algorithm does not converge (according to Heidelberger and Welch's convergence diagnostic), 
#' we automatically double N_MCMC and burn_in, and run it a second time (a message will be printed on screen to inform users).
#' @param burn_in the length of the burn-in; i.e., the initial part of the MCMC chain to be discarded (before convergence is reached).
#' Min 500.
#' If no convergence is reached, the 'burn_in' is authomatically increased (up to N_MCMC/2) according to 
#' the convergence detected by Heidelberger and Welch's convergence diagnostic.
#' If our algorithm does not converge even after increasing the burn-in, 
#' we automatically double N_MCMC and burn_in, and run it a second time (a message will be printed on screen to inform users).
#' 
#' @return A \code{list} of 3 \code{data.frame} objects.
#' 'Differential_results' contains results from differential testing only;
#' 'US_results' has results for the proportion of Spliced and Unspliced counts 
#' (Ambiguous counts are allocated 50:50 to Spliced and Unspliced);
#' 'USA_results' has results for the proportion of Spliced, Unspliced and Ambiguous counts 
#' (Ambiguous counts are reported separately from Spliced and Unspliced counts).
#' Columns 'Gene_id' and 'Cluster_id' contain the gene and cell-cluster name, 
#' while 'p_val', 'p_adj.loc' and 'p_adj.glb' report the raw p-values, locally and globally adjusted p-values, 
#' via Benjamini and Hochberg (BH) correction.
#' In locally adjusted p-values ('p_adj.loc') BH correction is applied to each cluster separately, 
#' while in globally adjusted p-values ('p_adj.glb') BH correction is performed to the results from all clusters.
#' Columns 'pi' and 'sd' indicate the proportion and standard deviation, respectively, 
#' 'S', 'U' and 'A' refer to Spliced, Unspliced and Ambiguous counts, respectively,
#' while 'gr_A' and 'gr_B' refer to group A and B, respectively.
#' For instance, columns 'pi_S-gr_A' and 'sd_S-gr_A' indicate the estimates and standard deviation (sd) 
#' for the proportion of Spliced (pi_S) and Unspliced (pi_U) counts in group A, respectively.
#' 
#' @examples
#' # load internal data to the package:
#' data_dir = system.file("extdata", package = "DifferentialRegulation")
#' 
#' # specify samples ids:
#' sample_ids = paste0("organoid", c(1:3, 16:18))
#' # set directories of each sample input data (obtained via alevin-fry):
#' base_dir = file.path(data_dir, "alevin-fry", sample_ids)
#' file.exists(base_dir)
#' 
#' # set paths to USA counts, cell id and gene id:
#' # Note that alevin-fry needs to be run with '--use-mtx' option
#' # to store counts in a 'quants_mat.mtx' file.
#' path_to_counts = file.path(base_dir,"/alevin/quants_mat.mtx")
#' path_to_cell_id = file.path(base_dir,"/alevin/quants_mat_rows.txt")
#' path_to_gene_id = file.path(base_dir,"/alevin/quants_mat_cols.txt")
#'
#' # load USA counts:
#' sce = load_USA(path_to_counts,
#'                path_to_cell_id,
#'                path_to_gene_id,
#'                sample_ids)
#'  
#' # define the design of the study:
#' design = data.frame(sample = sample_ids,
#'                     group = c( rep("3 mon", 3), rep("6 mon", 3) ))
#' design
#' 
#' # cell types should be assigned to each cell;
#' # here we load pre-computed cell types:
#' path_to_DF = file.path(data_dir,"DF_cell_types.txt")
#' DF_cell_types = read.csv(path_to_DF, sep = "\t", header = TRUE)
#' matches = match(colnames(sce), DF_cell_types$cell_id)
#' sce$cell_type = DF_cell_types$cell_type[matches]
#'
#' PB_counts = compute_PB_counts(sce = sce,
#'                               EC_list = NULL,
#'                               design =  design,
#'                               sample_col_name = "sample",
#'                               group_col_name = "group",
#'                               sce_cluster_name = "cell_type",
#'                               min_cells_per_cluster = 100, 
#'                               min_counts_per_gene_per_group = 20)
#'                               
#' # Differential regulation test based on estimated USA (unspliced, spliced, ambiguous) counts
#' set.seed(169612)
#' results_USA = DifferentialRegulation(PB_counts, EC = FALSE)
#' 
#' # DifferentialRegulation returns of a list of 3 data.frames:
#' # "Differential_results" contains results from differential testing only;
#' # "US_results" has estimates and standard deviation (SD) for pi_S and pi_U (proportion of Spliced and Unspliced counts);
#' # "USA_results" has estimates and standard deviation (SD) for pi_S, pi_U and pi_A (proportion of Spliced, Unspliced and Ambiguous counts).
#' names(results_USA)
#' 
#' # We visualize differential results:
#' head(results_USA$Differential_results)
#' 
#' # For improved performance, at a higher computational cost,
#' # we recommend using equivalence classes (EC) (here not run for computational reasons)
#' if(FALSE){
#'   # set paths to EC counts and ECs:
#'   path_to_EC_counts = file.path(base_dir,"/alevin/geqc_counts.mtx")
#'   path_to_EC = file.path(base_dir,"/alevin/gene_eqclass.txt.gz")
#' 
#'   # load EC counts:
#'   EC_list = load_EC(path_to_EC_counts,
#'                     path_to_EC,
#'                     path_to_cell_id,
#'                     path_to_gene_id,
#'                     sample_ids)
#'                     
#'   PB_counts = compute_PB_counts(sce = sce,
#'                                 EC_list = EC_list,
#'                                 design =  design,
#'                                 sample_col_name = "sample",
#'                                 group_col_name = "group",
#'                                 sce_cluster_name = "cell_type",
#'                                 min_cells_per_cluster = 100, 
#'                                 min_counts_per_gene_per_group = 20)
#'   
#'   # to reduce memory usage, we can remove the EC_list object:
#'   rm(EC_list)
#'   
#'   set.seed(169612) 
#'   results_EC = DifferentialRegulation(PB_counts)
#'   
#'   names(results_EC)
#'   
#'   # We visualize differential results:
#'   head(results_EC$Differential_results)
#' }
#' 
#' # plot top (i.e., most significant) result:
#' # plot USA proportions:
#' plot_pi(results_USA,
#'         type = "USA",
#'         gene_id = results_USA$Differential_results$Gene_id[1],
#'         cluster_id = results_USA$Differential_results$Cluster_id[1])
#' 
#' # plot US proportions:
#' plot_pi(results_USA,
#'         type = "US",
#'         gene_id = results_USA$Differential_results$Gene_id[1],
#'         cluster_id = results_USA$Differential_results$Cluster_id[1])
#'
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{load_EC}}, \code{\link{load_USA}}, \code{\link{plot_pi}}, 
#' 
#' @export
DifferentialRegulation = function(PB_counts,
                                  EC = TRUE,
                                  n_cores = NULL, # by default = n_clusters
                                  N_MCMC = 2000,
                                  burn_in = 500){
  if(N_MCMC < 2*10^3){
    message("'N_MCMC' must be at least 2*10^3")
    return(NULL)
  }
  
  if(burn_in < 500){
    message("'burn_in' must be at least 500")
    return(NULL)
  }
  
  # retrieve objects needed:
  PB_data_prepared = PB_counts[[1]]
  min_counts_per_gene_per_group = PB_counts[[2]]
  n_samples = PB_counts[[3]]
  n_samples_per_group = PB_counts[[4]]
  numeric_groups = PB_counts[[5]]
  cluster_ids_kept = PB_counts[[6]]
  sample_ids_per_group = PB_counts[[7]]
  n_groups = PB_counts[[8]]
  gene_ids_sce = PB_counts[[9]]
  n_cell_types = PB_counts[[10]]
  levels_groups = PB_counts[[11]]
  
  if(!EC){
    message("'EC' was set to 'FALSE': estimated counts will be used to perform differential testing (faster, but marginally less accurate).")
    message("We recommend using equivalence classes counts (slower, but marginally more accurate).")
  }else{
    min_counts_ECs = PB_counts[[12]]
    list_EC_gene_id = PB_counts[[13]]
    list_EC_USA_id = PB_counts[[14]]
    genes = PB_counts[[15]]
    n_genes = PB_counts[[16]]
  }
  
  rm(PB_counts)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # Register parallel cores:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(is.null(n_cores)){
    message("'n_cores' was left 'NULL'.")
    message("Since tasks are paralellized on cell clusters, we will set 'n_cores' to the number of clusters that will be analyzed.")
    message("'n_cores' set to: ", n_cell_types)
    n_cores = n_cell_types
  }
  
  # check if n_cores is the same length as n_cell_types:
  if(n_cores != n_cell_types){
    message("We detected ", n_cell_types, " cell clusters, while 'n_cores' was equal to ", n_cores)
    message("Since tasks are paralellized on cell clusters, we recommend setting 'n_cores' to the number of clusters")
    
    cores_equal_clusters = FALSE
  }else{
    cores_equal_clusters = TRUE
  }
  
  cl = makeCluster(n_cores, setup_strategy = "sequential")
  registerDoParallel(cl, n_cores)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # run MCMC in parallel:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(EC){
    RESULTS = MCMC_ECs(PB_data_prepared,
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
                       cores_equal_clusters)
  }else{
    RESULTS = MCMC_USA(PB_data_prepared,
                       min_counts_per_gene_per_group,
                       N_MCMC,
                       burn_in,
                       n_samples,
                       n_samples_per_group,
                       numeric_groups,
                       cl,
                       cluster_ids_kept,
                       sample_ids_per_group,
                       n_groups,
                       gene_ids_sce,
                       cores_equal_clusters)
  }
  
  #Error in { : task 5 failed - "object 'PB_data_prepared' not found"
    
  stopCluster(cl) 
  stopImplicitCluster()
  
  # separate
  RES = RESULTS[[1]]
  Convergence_results = RESULTS[[2]]
  rm(RESULTS)
  
  # CHECK if ALL NULL (if all clusters return null):
  if(is.null(RES)){
    message("Results are all NULL. This may be due to:")
    message("i) failure of convergence (check 'Convergence_results' object), or")
    message("ii) no gene-cluster combinations passing the minimum expression thresholds.")
    
    res = list( Differential_results = NULL,
                US_results = NULL,
                USA_results = NULL,
                Convergence_results = Convergence_results)
    
    return(res)
  }
  
  # Replace group names, with those provided in design
  names = colnames(RES)
  sel_A = grep("gr_A", names )
  sel_B = grep("gr_B", names )
  
  colnames(RES)[sel_A] = gsub("gr_A", levels_groups[1], names[sel_A] )
  colnames(RES)[sel_B] = gsub("gr_B", levels_groups[2], names[sel_B] )
  
  # order results by significance (raw p-value)
  ord = order(RES$p_val)
  if(! ( any(is.na(ord)) | any(is.null(ord)) | any(is.nan(ord)) ) ){
    RES = RES[ ord, ]
  }
  
  res = list( Differential_results = RES[, seq_len(5)],
              US_results = RES[, seq_len(13)],
              USA_results = RES[, c(seq.int(1,5,by = 1),seq.int(14,25,by = 1))],
              Convergence_results = Convergence_results)
  
  res
}