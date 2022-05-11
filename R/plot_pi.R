#' Plot the estimated proportions of US or USA counts in each group
#'
#' \code{plot_pi} plots the posterior means of the proportions
#' of US (if 'type' = 'US') or USA (if 'type' = 'USA') counts, in each group.
#' If 'CI' is TRUE, a profile Wald type confidence interval will also be added;
#' the level of the confidence interval is specified by 'CI_level'.
#' 
#' @param results a \code{list} of 3 \code{\linkS4class{data.frame}} objects, computed via \code{DifferentialRegulation}.
#' @param gene_id a character, indicating the gene to plot.
#' @param cluster_id a character, indicating the cell cluster to plot.
#' @param type a character (either 'SU' or 'SUA').
#' If "SU", it plots the proportion of Spliced and Unspliced counts 
#' (Ambiguous counts are assigned 50:50 to Spliced and Unspliced counts).
#' If "SUA" (default), it plots the proportion of Spliced, Unspliced and Ambiguous counts 
#' (Ambiguous counts are kept separately).
#' Note that, although US reads are easier to interpret, 
#' USA reads more closely reflect what is being compared between conditions.
#' @param CI a logical ('TRUE' by default), indicating whether to plot a 
#' profile Wald type confidence interval around the estimated proportions.
#' @param CI_level a numeric between 0 and 1, indicating the level of the confidence interval.
#' @return A \code{ggplot} object.
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
#' # DifferentialRegulation returns of a list of 3 data.frames:
#' # "Differential_results" contains results from differential testing only;
#' # "US_results" has estimates and standard deviation (SD) for pi_S and pi_U (proportion of Spliced and Unspliced counts);
#' # "USA_results" has estimates and standard deviation (SD) for pi_S, pi_U and pi_A (proportion of Spliced, Unspliced and Ambiguous counts).
#' names(results_USA)
#' 
#' head(results_USA[[1]])
#' head(results_USA[[2]])
#' 
#' # We can also sort results by significance, if we want, before visualizing them.
#' DR = results_USA[[1]]
#' DR = DR[ order(DR$p_val), ]
#' head(DR)
#' 
#' # For improved performance, at a higher computational cost,
#' # we recommend using equivalence classes (EC) (here not run for computational reasons)
#' # see help(DifferentialRegulation) examples.
#' 
#' # plot top (i.e., most significant) result:
#' # plot USA proportions:
#' plot_pi(results_USA,
#'         type = "USA",
#'         gene_id = DR$Gene_id[1],
#'         cluster_id = DR$Cluster_id[1])
#' 
#' # plot US proportions:
#' plot_pi(results_USA,
#'         type = "US",
#'         gene_id = DR$Gene_id[1],
#'         cluster_id = DR$Cluster_id[1])
#'
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{DifferentialRegulation}}
#' 
#' @export
plot_pi = function(results,
                   gene_id,
                   cluster_id,
                   type = "USA",
                   CI = TRUE,
                   CI_level = 0.95){
  if( (!is.list(results)) | (length(results) != 3) ){
    message("'results' must be a 'list' of 3 objects, as created via 'DifferentialRegulation' function.")
    return(NULL)
  }
  if(!all(c("US_results", "USA_results") %in% names(results))){
    message("'names(results)' must include both '' and 'USA_results'.")
    message("'results' must be a 'list' of 3 objects, as created via 'DifferentialRegulation' function.")
    return(NULL)
  }
  
  if(!is.logical(CI)){
    message("'CI' must be 'TRUE' or 'FALSE'")
    return(NULL)
  }
  
  if(!is.numeric(CI_level)){
    message("'CI_level' must be 'numeric'.")
    return(NULL)
  }
  if( (CI_level < 0) | (CI_level > 1) ){
    message("'CI_level' must be between 0 and 1.")
    return(NULL)
  }
  
  if(!is.character(type)){
    message("'type' must be a 'character'.")
    return(NULL)
  }
  if( length(type) != 1){
    message("'type' contains ", length(type), " values; it should only include 1: either 'US' or 'USA'.")
    return(NULL)
  }
  if( ! type %in% c("US", "USA") ){
    message("'type' is ", type)
    message("'type' must be use either 'US' or 'USA'.")
    return(NULL)
  }
  
  if(!is.character(gene_id)){
    gene_id = as.character(gene_id)
  }
  if( length(gene_id) != 1){
    message("'gene_id' contains ", length(gene_id), " values: one gene only can be specified.")
    return(NULL)
  }
  
  if(!is.character(cluster_id)){
    gene_id = as.character(gene_id)
  }
  if( length(cluster_id) != 1){
    message("'cluster_id' contains ", length(cluster_id), " values: one cluster_id only can be specified.")
    return(NULL)
  }
  
  #select DF and cols for pi and SD
  if(type == "US"){
    DF = results$US_results
    sel_pi = 6:9
    sel_sd = 10:13
    tr_names = factor(c("S", "U"), levels = c("S", "U") )
  }else{
    DF = results$USA_results
    sel_pi = 6:11
    sel_sd = 12:17
    tr_names = factor(c("S", "U", "A"), levels = c("S", "U", "A") )
  }
  
  sel = which( (DF$Gene_id == gene_id) & (DF$Cluster_id == cluster_id) )
  
  if( length(sel) == 0 ){
    if(type == "US"){
      message("'gene_id' and 'cluster_id' combination not found in 'results$US_results$Gene_id'")
    }else{
      message("'gene_id' and 'cluster_id' combination not found in 'results$USA_results$Gene_id'")
    }
    return(NULL)
  }
  
  #group_names = levels(factor(x@samples_design$group))
  group_names = c("A", "B")
  
  pi = DF[sel, sel_pi]
  SD =  DF[sel, sel_sd]
  n_groups = 2
  
  prop_samp = data.frame(feature_id = factor( rep(tr_names, n_groups)), 
                         proportion = unlist(c(pi)),
                         LB = pmax(0, unlist(c(pi - qnorm(1 - (1-CI_level)/2) * SD)) ), # LB must be >= 0
                         UB = pmin(1, unlist(c(pi + qnorm(1 - (1-CI_level)/2) * SD)) ), # UB must be <= 0
                         group = rep(group_names, each = length(tr_names)),
                         stringsAsFactors = FALSE)
  
  # Plot the estimated average proportions of each groups:
  ggp = ggplot() +
    geom_bar(data = prop_samp, aes_string(x = "feature_id", y = "proportion", 
                                          fill = "group"),
             stat = "identity", position = position_dodge(width = 0.9)) +
    theme_bw() + 
    theme(axis.text.x = element_text(vjust = 0.5), 
          axis.text = element_text(size=16), 
          axis.title = element_text(size=14, face="bold"), 
          plot.title = element_text(size=16), 
          legend.position = "right", 
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14)) +
    ggtitle(paste("Gene:", gene_id, "- Cluster:", cluster_id)) +
    xlab("Features") +
    ylab("Proportions") +
    geom_point(data = prop_samp, 
               aes_string(x = "feature_id", y = "proportion", 
                          group = "group", fill = "group"), 
               position = position_dodge(width = 0.9), size = 3, shape = 1, 
               alpha = 0.75)
  
  if(CI){
    ggp = ggp +
      geom_errorbar(data = prop_samp,
                    aes_string(x = "feature_id", ymin = "LB", ymax = "UB",
                               group = "group"), 
                    position = position_dodge(width = 0.9), size = 0.5, 
                    width = 0.5,
                    alpha = 0.5)
  }
  
  ggp
}