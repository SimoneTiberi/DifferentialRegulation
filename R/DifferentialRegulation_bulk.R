#' Discover differentially regulated genes from bulk RNA-seq data
#'
#' \code{DifferentialRegulation_bulk} identified differentially regulated genes between two conditions 
#' (e.g., healthy vs. disease or treated vs. untreated) in each cluster of cells.
#' Parameters are inferred via Markov chain Monte Carlo (MCMC) techniques and a differential testing is performed 
#' via a multivariate Wald test on the posterior densities of the group-level US (Unspliced, and Spliced) 
#' counts relative abundance.
#' 
#' @param SE a \code{SummarizedExperiment} object, computed via \code{\link{load_bulk_US}}.
#' @param EC_list a \code{list}, computed via \code{\link{load_bulk_EC}}.
#' @param design a \code{\linkS4class{data.frame}} indicating the design of the experiment with one row for each sample;
#' 'design' must contain a column with the sample id and one with the group id.
#' @param sample_col_name a character ("sample" by default), indicating the column name of the 'design' element which stores the sample id.
#' @param group_col_name a character ("group" by default), indicating the column name of the 'design' element which stores the group id.
#' @param min_counts_per_transcript_per_group minimum number of counts per transcript, across all samples in a group
#' Only genes with at least 'min_counts_per_transcript_per_group' counts in both groups of samples will be analyzed.
#' @param min_counts_ECs equivalence classes (ECs) filter.
#' 'min_counts_ECs' indicates the minimum number of counts (across all cells in a cell cluster) for each equivalence class;
#' by default all ECs are considered (min_counts_ECs = 0).
#' ECs with less or equal than 'min_counts_ECs' will be discarded.
#' Increasing 'min_counts_ECs' will marginally decrease computational cost computational at the cost of a marginal loss in performance.
#' @param n_cores the number of cores to parallelize the tasks on.
#' Note that only a minor part of the function runs in parallel (i.e., the prior for the dispersion),
#' while most of the function does not (e.g., the MCMC).
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
#' @param c_prop temporary parameter; do NOT edit.
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
#' sample_ids = paste0("sample", seq_len(6))
#' 
#' # US estimates:
#' quant_files = file.path(data_dir, "salmon", sample_ids, "quant.sf")
#' file.exists(quant_files)
#' 
#' # Equivalence classes:
#' equiv_classes_files = file.path(data_dir, "salmon", sample_ids, "aux_info/eq_classes.txt.gz")
#' file.exists(equiv_classes_files)
#' 
#' # load EC:
#' EC_list = load_bulk_EC(path_to_eq_classes = equiv_classes_files,
#'                        n_cores = 2)
#' 
#' # load US estimated counts:
#' SE = load_bulk_US(quant_files,
#'                   sample_ids)
#' 
#' # define the design of the study:
#' group_names = rep(c("A", "B"), each = 3)
#' design = data.frame(sample = sample_ids,
#'                     group = group_names)
#' design
#' 
#' set.seed(169612) 
#' results = DifferentialRegulation_bulk(SE = SE, 
#'                                       EC_list = EC_list,
#'                                       design = design, 
#'                                       n_cores = 2)
#' 
#' names(results)
#'   
#' # We visualize differential results:
#' head(results$Differential_results)
#' 
#' # plot top (i.e., most significant) result:
#' # plot USA proportions:
#' plot_bulk_pi(results,
#'              transcript_id = results$Differential_results$Transcript_id[1])
#' 
#' @author Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{load_bulk_EC}}, \code{\link{load_bulk_US}}, \code{\link{plot_pi}}
#' 
#' @export
DifferentialRegulation_bulk = function(SE,
                                       EC_list,
                                       design,
                                       sample_col_name = "sample",
                                       group_col_name = "group",
                                       min_counts_per_transcript_per_group = 10,
                                       N_MCMC = 2000,
                                       burn_in = 500,
                                       min_counts_ECs = 0,
                                       n_cores = 1,
                                       undersampling_int = 10,
                                       c_prop = 0.3){
  if( (undersampling_int < 1) | (undersampling_int > 10) ){
    message("'undersampling_int' must be an integer between 1 and 10 (included)")
    return(NULL)
  }
  if( round(undersampling_int) != undersampling_int ){
    message("'undersampling_int' must be an integer")
    return(NULL)
  }
  
  if(!is.numeric(min_counts_per_transcript_per_group)){
    message("'min_counts_per_transcript_per_group' must be numeric.")
  }
  if(min_counts_per_transcript_per_group < 10){
    message("'min_counts_per_transcript_per_group' must be at least 10.")
  }
  if(!is.numeric(min_counts_ECs)){
    message("'min_counts_ECs' must be numeric.")
  }
  
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
  
  groups = factor(design[, sel_col ])
  levels_groups = levels(groups)
  n_groups = length(levels_groups)
  numeric_groups = as.numeric(groups)
  
  if(n_groups > 2.5){
    message("We detected ", n_groups, " groups in the design")
    message("At present, only 2 group comparisons are implemented")
    return(NULL)
  }
  
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
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # get pseudo-bulk EC counts and USA counts from SE:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # filter lowly abundant transcripts
  tot = assays(SE)$spliced + assays(SE)$unspliced
  sel = ( rowSums(tot[ ,groups == levels_groups[1]]) >= min_counts_per_transcript_per_group ) & 
    ( rowSums(tot[ ,groups == levels_groups[2]]) >= min_counts_per_transcript_per_group )
  
  # filter transcripts with S version only:
  sel = sel & rowData(SE)$has_U_tr
  
  message(sum(sel)," transcripts pass filtering and will be analyzed.")
  if( sum(sel) == 0){
    return(NULL)
  }
  
  # get all transcript names:
  tr_ids_sce = rownames(SE)
  # get transcript names that pass filters
  tr_ids_sce_pass_filter = tr_ids_sce[sel]
  
  res = MCMC_bulk_EC(SE,
                     EC_list,
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
                     c_prop)
  
  res
}