#' Plot the estimated proportions of US or USA counts in each group - single-cell RNA-seq data
#'
#' \code{plot_pi} plots the posterior means of the proportions
#' of US (if 'type' = 'US') or USA (if 'type' = 'USA') counts, in each group.
#' If 'CI' is TRUE, a profile Wald type confidence interval will also be added;
#' the level of the confidence interval is specified by 'CI_level'.
#' 
#' @param results a \code{list} of \code{\linkS4class{data.frame}} objects, computed via \code{\link{DifferentialRegulation}}.
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
#' # see the example of DifferentialRegulation function:
#' help(DifferentialRegulation)
#' 
#' @author Simone Tiberi \email{simone.tiberi@unibo.it}
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
  if( !is.list(results) ){
    message("'results' must be a 'list' created via 'DifferentialRegulation' function.")
    return(NULL)
  }
  if(!all(c("US_results", "USA_results") %in% names(results))){
    message("'names(results)' must include both '' and 'USA_results'.")
    message("'results' must be a 'list' created via 'DifferentialRegulation' function.")
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
    sel_pi = 7:10
    sel_sd = 11:14
    tr_names = factor(c("S", "U"), levels = c("S", "U") )
    
    group_names = substring(colnames( DF )[c(7,9)],6)
  }else{
    DF = results$USA_results
    sel_pi = 7:12
    sel_sd = 13:18
    tr_names = factor(c("S", "U", "A"), levels = c("S", "U", "A") )
    
    group_names = substring(colnames( DF )[c(7,10)],6)
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
                    position = position_dodge(width = 0.9), linewidth = 0.5, 
                    width = 0.5,
                    alpha = 0.5)
  }
  
  ggp
}