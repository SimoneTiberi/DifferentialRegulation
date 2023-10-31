#' Discover differentially regulated genes from single-cell RNA-seq data
#'
#' \code{DifferentialRegulation} identified differentially regulated genes between two conditions 
#' (e.g., healthy vs. disease or treated vs. untreated) in each cluster of cells.
#' Parameters are inferred via Markov chain Monte Carlo (MCMC) techniques and a differential testing is performed 
#' via a multivariate Wald test on the posterior densities of the group-level USA (Unspliced, Spliced and Ambiguous) counts relative abundance.
#' 
#' @param PB_counts a \code{list}, computed via \code{\link{compute_PB_counts}}
#' @param n_cores the number of cores to parallelize the tasks on.
#' Since parallelization is at the cluster level (each cluster is parallelized on a thread), 
#' we suggest setting n_cores to the number of clusters (e.g., cell-types), as set by default if 'n_cores' is not specified.
#' @param N_MCMC the number of iterations for the MCMC algorithm (including burn-in). Min 2*10^3.
#' If our algorithm does not converge (according to Heidelberger and Welch's convergence diagnostic), 
#' we automatically double N_MCMC and burn_in, and run it a second time (a message will be printed on screen to inform users).
#' @param burn_in the length of the burn-in; i.e., the initial part of the MCMC chain to be discarded (before convergence is reached).
#' Min 500.
#' If no convergence is reached, the 'burn_in' is automatically increased (up to N_MCMC/2) according to 
#' the convergence detected by Heidelberger and Welch's convergence diagnostic.
#' If our algorithm does not converge even after increasing the burn-in, 
#' we automatically double N_MCMC and burn_in, and run it a second time (a message will be printed on screen to inform users).
#' @param undersampling_int the undersampling of the latent variables.
#' While model parameters are sampled at each iteration,
#' RNA-seq counts are allocated to their transcript (and spliced/unspliced) version of origin,
#' every `undersampling_int` iterations.
#' Increasing `undersampling_int` will decrease the runtime, but may marginally affect performance.
#' In our benchmarks, no differences in performance were observed for values up to 10.
#' @param traceplot a logical value indicating whether to return the posterior chain
#' of "pi_U", for both groups (i.e., the group-level relative abundance of unspliced reads).
#' If TRUE, the posterior chains are stored in 'MCMC_U' object,
#' and can be plotted via 'plot_traceplot' function.
#' 
#' @return A \code{list} of 4 \code{data.frame} objects.
#' 'Differential_results' contains results from differential testing only;
#' 'US_results' has results for the proportion of Spliced and Unspliced counts 
#' (Ambiguous counts are allocated 50:50 to Spliced and Unspliced);
#' 'USA_results' includes results for the proportion of Spliced, Unspliced and Ambiguous counts 
#' (Ambiguous counts are reported separately from Spliced and Unspliced counts);
#' 'Convergence_results' contains information about convergence of posterior chains.
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
#' # set paths to EC counts and ECs:
#' path_to_EC_counts = file.path(base_dir,"/alevin/geqc_counts.mtx")
#' path_to_EC = file.path(base_dir,"/alevin/gene_eqclass.txt.gz")
#' 
#' # load EC counts:
#' EC_list = load_EC(path_to_EC_counts,
#'                   path_to_EC,
#'                   path_to_cell_id,
#'                   path_to_gene_id,
#'                   sample_ids)
#'                     
#' PB_counts = compute_PB_counts(sce = sce,
#'                               EC_list = EC_list,
#'                               design =  design,
#'                               sample_col_name = "sample",
#'                               group_col_name = "group",
#'                               sce_cluster_name = "cell_type",
#'                               min_cells_per_cluster = 100, 
#'                               min_counts_per_gene_per_group = 20)
#' 
#' # to reduce memory usage, we can remove the EC_list object:
#' rm(EC_list)
#'   
#' set.seed(1609612) 
#' results = DifferentialRegulation(PB_counts,
#'                                  n_cores = 2,
#'                                  traceplot = TRUE)
#'   
#' names(results)
#'   
#' # We visualize differential results:
#' head(results$Differential_results)
#' 
#' # plot top (i.e., most significant) result:
#' # plot USA proportions:
#' plot_pi(results,
#'         type = "USA",
#'         gene_id = results$Differential_results$Gene_id[1],
#'         cluster_id = results$Differential_results$Cluster_id[1])
#' 
#' # plot US proportions:
#' plot_pi(results,
#'         type = "US",
#'         gene_id = results$Differential_results$Gene_id[1],
#'         cluster_id = results$Differential_results$Cluster_id[1])
#'        
#' # plot the corresponding traceplot:
#' plot_traceplot(results,
#'                gene_id = results$Differential_results$Gene_id[1],
#'                cluster_id = results$Differential_results$Cluster_id[1])
#'
#' @author Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{load_EC}}, \code{\link{load_USA}}, code{\link{compute_PB_counts}}, \code{\link{plot_pi}}, \code{\link{plot_traceplot}}
#' 
#' @export
DifferentialRegulation = function(PB_counts,
                                  n_cores = NULL, # by default = n_clusters
                                  N_MCMC = 2000,
                                  burn_in = 500,
                                  undersampling_int = 10,
                                  traceplot = FALSE){
  if( (undersampling_int < 1) | (undersampling_int > 10) ){
    message("'undersampling_int' must be an integer between 1 and 10 (included)")
    return(NULL)
  }
  if( round(undersampling_int) != undersampling_int ){
    message("'undersampling_int' must be an integer")
    return(NULL)
  }
  
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
  
  min_counts_ECs = PB_counts[[12]]
  list_EC_gene_id_original = PB_counts[[13]]
  list_EC_USA_id_original = PB_counts[[14]]
  genes = PB_counts[[15]]
  n_genes = PB_counts[[16]]
  
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
  if(n_cores < n_cell_types){
    message("We detected ", n_cell_types, " cell clusters, while 'n_cores' was equal to ", n_cores)
    message("Since tasks are paralellized on cell clusters, we recommend setting 'n_cores' to the number of clusters")
    
    cores_equal_clusters = FALSE
  }else{
    cores_equal_clusters = TRUE
  }
  
  # if at least 2 clusters:
  cluster = makeCluster(n_cores, setup_strategy = "sequential")
  registerDoParallel(cluster, n_cores)
  # pass libPath to workers:
  clusterCall(cluster, function(x) .libPaths(x), .libPaths())
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # run MCMC in parallel:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  RESULTS = MCMC_ECs(PB_data_prepared,
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
                     traceplot)
  
  stopCluster(cluster) 
  stopImplicitCluster()
  
  # separate
  RES = RESULTS[[1]]
  Convergence_results = RESULTS[[2]]
  MCMC_U = RESULTS[[3]]
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
  
  if(traceplot){
    res = list( Differential_results = RES[, seq_len(6)],
                US_results = RES[, seq_len(14)],
                USA_results = RES[, c(seq.int(1,6,by = 1),seq.int(15,26,by = 1))],
                Convergence_results = Convergence_results,
                MCMC_U = MCMC_U)
  }else{
    res = list( Differential_results = RES[, seq_len(6)],
                US_results = RES[, seq_len(14)],
                USA_results = RES[, c(seq.int(1,6,by = 1),seq.int(15,26,by = 1))],
                Convergence_results = Convergence_results)
    
  }
  
  res
}