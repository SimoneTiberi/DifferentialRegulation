#' Discover differentially regulated genes
#'
#' \code{DifferentialRegulation} identified differentially regulated genes between two conditions 
#' (e.g., healthy vs. disease or treated vs. untreated) in each cluster of cells.
#' Parameters are inferred via Markov chain Monte Carlo (MCMC) techniques and a differential testing is performed 
#' via a multivariate Wald test on the posterior densities of the group-level USA (Unspliced, Spliced and Ambiguous) counts relative abundance.
#' 
#' @param sce a \code{SingleCellExperiment} object, computed via \code{load_USA}.
#' @param EC_list a \code{list}, computed via \code{load_EC}.
#' @param samples_design TODO
#' @param sample_col_name TODO
#' @param group_col_name TODO
#' @param sce_cluster_name TODO
#' @param min_cells_per_cluster TODO
#' @param min_counts_per_gene_per_group TODO
#' @param n_cores TODO
#' @param N_MCMC TODO
#' @param burn_in TODO
#' @param min_counts_ECs TODO
#' 
#' @return A \code{data.frame} object.
#' 
#' @examples
#' TODO
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{load_EC}}, \code{\link{load_USA}}
#' 
#' @export
DifferentialRegulation = function(sce,
                                  EC_list = NULL,
                                  samples_design,
                                  sample_col_name = "sample",
                                  group_col_name = "group",
                                  sce_cluster_name = "cell_type",
                                  min_cells_per_cluster = 100, 
                                  min_counts_per_gene_per_group = 20,
                                  n_cores = NULL, # by default = n_clusters
                                  N_MCMC = 2000,
                                  burn_in = 500,
                                  min_counts_ECs = 0
){
  if( !is.data.frame(samples_design) ){
    message("'samples_design' must be a data.frame object")
    return(NULL)
  }
  # select the column of samples_design which is called 'group_col_name'
  if( !(group_col_name %in% colnames(samples_design)) ){
    message("Column ", group_col_name, " missing in 'samples_design'")
    message("'group_col_name' should specify the column name of 'samples_design' containing the group id of each sample")
    return(NULL)
  }
  sel_col = which(group_col_name == colnames(samples_design))
  if( length(sel_col) > 1.5 ){
    message( length(sel_col) , " columns from 'samples_design' are called ", group_col_name)
    message("Remove duplicated columns from 'samples_design' and provide a unique column for the group id")
    return(NULL)
  }
  groups = factor(samples_design[, sel_col ])
  levels_groups = levels(groups)
  n_groups = length(levels_groups)
  numeric_groups = as.numeric(groups)
  
  sample_ids_per_group = lapply(1:n_groups, function(gr){
    which(numeric_groups == gr) - 1 # -1 !
  })
  n_samples_per_group = sapply(sample_ids_per_group, length)
  
  # select the column of samples_design which is called 'sample_col_name'
  if( !(sample_col_name %in% colnames(samples_design)) ){
    message("Column ", sample_col_name, " missing in 'samples_design'")
    message("'sample_col_name' should specify the column name of 'samples_design' containing the group id of each sample")
    return(NULL)
  }
  sel_col = which(sample_col_name == colnames(samples_design))
  if( length(sel_col) > 1.5 ){
    message( length(sel_col) , " columns from 'samples_design' are called ", sample_col_name)
    message("Remove duplicated columns from 'samples_design' and provide a unique column for the group id")
    return(NULL)
  }
  samples = samples_design[, sel_col ]
  n_samples = length(samples)
  
  # cluster ids:
  sel = which(names(colData(x)) == sce_cluster_name)
  if( length(sel) == 0 ){
    message("'sce_cluster_name' not found in names(colData(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'sce_cluster_name' found multiple times in names(colData(x))")
    return(NULL)
  }
  clusters = factor(colData(x)[[sel]])
  n_clusters = nlevels(clusters)
  # clusters = as.integer(as.numeric(clusters)-1)
  
  # select cell types with at least xx cells across all samples
  table_clusters = table(clusters)
  cluster_ids_kept = names(table_clusters[table_clusters >= min_cells_per_cluster])
  
  message("the following cell types have more than", min_cells_per_cluster, "cells and will be analyzed:")
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
  list_EC_gene_id = EC_list[[2]]
  list_EC_USA_id = EC_list[[3]]
  genes = EC_list[[4]]
  n_genes = length(genes)
  
  gene_ids_sce = rownames(sce)
  
  # rm (heavy) unnecessary objects:
  rm(sce); rm(clusters); rm(EC_list)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # Register parallel cores:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(is.null(n_cores)){
    n_cores = n_cell_types
  }
  
  suppressWarnings({
    cl = makeCluster(n_cores, setup_strategy = "sequential")
  })
  registerDoParallel(cl, n_cores);
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # run MCMC in parallel:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  message("Starting the MCMC")
  
  if(is.null(EC_list)){
    message("'EC_list' was not provided: estimated counts in 'sce' will be used to perform differential testing (faster, but marginally less accurate).")
    message("If you want to use equivalence classes counts (recommended option: slower, but marginally more accurate), provide an 'EC_list' object, computed via 'load_EC' function.")
    
    RES = MCMC_sce(PB_data_prepared,
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
                   gene_ids_sce)
  }else{
    RES = MCMC_ECs(PB_data_prepared,
                   samples_design,
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
                   gene_ids_sce)
  }
  stopCluster(cl) 
  stopImplicitCluster()
  
  RES
}