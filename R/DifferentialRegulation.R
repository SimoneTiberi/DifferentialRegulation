#' Discover differentially regulated genes
#'
#' \code{DifferentialRegulation} identified differentially regulated genes between two conditions 
#' (e.g., healthy vs. disease or treated vs. untreated) in each cluster of cells.
#' Parameters are inferred via Markov chain Monte Carlo (MCMC) techniques and a differential testing is performed 
#' via a multivariate Wald test on the posterior densities of the group-level USA (Unspliced, Spliced and Ambiguous) counts relative abundance.
#' 
#' @param sce a \code{SingleCellExperiment} object, computed via \code{load_USA}.
#' @param EC_list a \code{list}, computed via \code{load_EC}.
#' @param design a \code{\linkS4class{data.frame}} indicating the design of the experiment with one row for each sample;
#' 'design' must contain a column with the sample id and one with the group id.
#' @param sample_col_name a character ("sample" by default), indicating the column name of the 'design' element which stores the sample id.
#' @param group_col_name a character ("group" by default), indicating the column name of the 'design' element which stores the group id.
#' @param sce_cluster_name a character ("cell_type" by default), indicating the name of the 'colData(sce)' element, 
#' which stores the cluster id of each cell (i.e., colData(sce)$name_cluster).
#' @param min_cells_per_cluster cell cluster (e.g., cell-type) filter.
#' 'min_cells_per_cluster' is the minimum number of cells, across all samples and groups, for a cell cluster to be considered.
#' Cell clusters with less than 'min_cells_per_cluster' cells will not be analyzed.  
#' @param min_counts_per_gene_per_group minimum number of counts per gene, in each cell, across all samples of every group.
#' In each cell cluster, only genes with at least 'min_counts_per_gene_per_group' counts in both groups of samples will be analyzed.  
#' @param min_counts_ECs equivalence classes (ECs) filter (NB: only used when 'EC_list' is provided)
#' 'min_counts_ECs' indicates the minimum number of counts (across all cells in a cell cluster) for each equivalence class;
#' by default all ECs are considered (min_counts_ECs = 0).
#' ECs with less or equal than 'min_counts_ECs' will be discarded.
#' Increasing 'min_counts_ECs' will marginally decrease computational cost computational at the cost of a marginal loss in performance.
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
#' @return A \code{data.frame} object.
#' 
#' @examples
#' # load internal data to the package:
#' data_dir = system.file("extdata", package = "DifferentialRegulation")
#' 
#' # specify 4 samples ids:
#' sample_ids = paste0("sample_", seq_len(4))
#' # set directories of each sample input data (obtained via alevin-fry):
#' base_dir = file.path(data_dir, "alevin-fry", sample_ids)
#' file.exists(base_dir)
#' 
#' # set paths to USA counts, cell id and gene id:
#' # Note that alevin-fry needs to be run with `--use-mtx` option
#' # to store counts in a `quants_mat.mtx` file.
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
#' group = c("A", "A", "B", "B"))
#' design
#' 
#' # cell types should be assigned to each cell;
#' # here we assume all cells belong to cell-type "cell":
#' sce$cell_type = "cell"
#'                
#' # Differential regulation test based on estimated USA (unspliced, spliced, ambiguous) counts
#' set.seed(169612)
#' results_USA = DifferentialRegulation(sce = sce,
#'                                     EC_list = NULL,
#'                                     design =  design,
#'                                     sample_col_name = "sample",
#'                                     group_col_name = "group",
#'                                     sce_cluster_name = "cell_type",
#'                                     min_cells_per_cluster = 100, 
#'                                     min_counts_per_gene_per_group = 20)
#' head(results_USA)
#' 
#' # We can also sort results by significance, if we want, before visualizing them.
#' results_USA = results_USA[ order(results_USA$p_val), ]
#' head(results_USA)
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
#'   results_EC = DifferentialRegulation(sce = sce,
#'                                       EC_list = EC_list,
#'                                       design =  design,
#'                                       sample_col_name = "sample",
#'                                       group_col_name = "group",
#'                                       sce_cluster_name = "cell_type",
#'                                       min_cells_per_cluster = 100, 
#'                                       min_counts_per_gene_per_group = 20)
#'   head(results_EC)
#' }
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{load_EC}}, \code{\link{load_USA}}
#' 
#' @export
DifferentialRegulation = function(sce,
                                  EC_list = NULL,
                                  design,
                                  sample_col_name = "sample",
                                  group_col_name = "group",
                                  sce_cluster_name = "cell_type",
                                  min_cells_per_cluster = 100, 
                                  min_counts_per_gene_per_group = 20,
                                  min_counts_ECs = 0,
                                  n_cores = NULL, # by default = n_clusters
                                  N_MCMC = 2000,
                                  burn_in = 500){
  if( !is.data.frame(design) ){
    message("'design' must be a data.frame object")
    return(NULL)
  }
  # select the column of design which is called 'group_col_name'
  if( !(group_col_name %in% colnames(design)) ){
    message("Column ", group_col_name, " missing in 'design'")
    message("'group_col_name' should specify the column name of 'design' containing the group id of each sample")
    return(NULL)
  }
  sel_col = which(group_col_name == colnames(design))
  if( length(sel_col) > 1.5 ){
    message( length(sel_col) , " columns from 'design' are called ", group_col_name)
    message("Remove duplicated columns from 'design' and provide a unique column for the group id")
    return(NULL)
  }
  
  if(N_MCMC < 2*10^3){
    message("'N_MCMC' must be at least 2*10^3")
    #return(NULL)
  }
  
  if(burn_in < 500){
    message("'burn_in' must be at least 500")
    #return(NULL)
  }
  
  groups = factor(design[, sel_col ])
  levels_groups = levels(groups)
  n_groups = length(levels_groups)
  numeric_groups = as.numeric(groups)
  
  sample_ids_per_group = lapply(seq_len(n_groups), function(gr){
    which(numeric_groups == gr) - 1 # -1 !
  })
  n_samples_per_group = vapply(sample_ids_per_group, length, FUN.VALUE = integer(1) )
  
  # select the column of design which is called 'sample_col_name'
  if( !(sample_col_name %in% colnames(design)) ){
    message("Column ", sample_col_name, " missing in 'design'")
    message("'sample_col_name' should specify the column name of 'design' containing the group id of each sample")
    return(NULL)
  }
  sel_col = which(sample_col_name == colnames(design))
  if( length(sel_col) > 1.5 ){
    message( length(sel_col) , " columns from 'design' are called ", sample_col_name)
    message("Remove duplicated columns from 'design' and provide a unique column for the group id")
    return(NULL)
  }
  samples = design[, sel_col ]
  n_samples = length(samples)
  
  # cluster ids:
  sel = which(names(colData(sce)) == sce_cluster_name)
  if( length(sel) == 0 ){
    message("'sce_cluster_name' not found in names(colData(sce))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'sce_cluster_name' found multiple times in names(colData(sce))")
    return(NULL)
  }
  clusters = factor(colData(sce)[[sel]])
  n_clusters = nlevels(clusters)
  # clusters = as.integer(as.numeric(clusters)-1)
  
  # select cell types with at least xx cells across all samples
  table_clusters = table(clusters)
  cluster_ids_kept = names(table_clusters[table_clusters >= min_cells_per_cluster])
  
  message("the following cell types have more than ", min_cells_per_cluster, " cells and will be analyzed:")
  message(cluster_ids_kept)
  
  n_cell_types = length(cluster_ids_kept)
  if( n_cell_types == 0 ){
    return(NULL)
  }
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # get pseudo-bulk EC counts and USA counts from sce:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  PB_data_prepared = lapply(cluster_ids_kept,
                            prepare_PB_counts,
                            sce = sce, clusters = clusters,
                            n_samples = n_samples, EC_list = EC_list)
  
  gene_ids_sce = rownames(sce)
  
  # rm (heavy) unnecessary objects:
  rm(sce); rm(clusters);
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # Register parallel cores:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(is.null(n_cores)){
    n_cores = n_cell_types
  }
  if(n_cores != n_cell_types){
    message("We detected ", n_cell_types, " cell clusters, while 'n_cores' was equal to, ", n_cores)
    message("Since tasks are paralellized on cell clusters, we recommend setting 'n_cores' to the number of clusters")
  }
  
  cl = makeCluster(n_cores, setup_strategy = "sequential")
  registerDoParallel(cl, n_cores);
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # run MCMC in parallel:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(is.null(EC_list)){
    message("'EC_list' was not provided: estimated counts in 'sce' will be used to perform differential testing (faster, but marginally less accurate).")
    message("If you want to use equivalence classes counts (recommended option: slower, but marginally more accurate), provide an 'EC_list' object, computed via 'load_EC' function.")
    
    RES = MCMC_sce(PB_data_prepared,
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
                   gene_ids_sce)
  }else{
    list_EC_gene_id = EC_list[[2]]
    list_EC_USA_id = EC_list[[3]]
    genes = EC_list[[4]]
    n_genes = length(genes)
    
    rm(EC_list) # rm after check!
    
    RES = MCMC_ECs(PB_data_prepared,
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
                   list_EC_USA_id)
  }
  stopCluster(cl) 
  stopImplicitCluster()
  
  RES
}