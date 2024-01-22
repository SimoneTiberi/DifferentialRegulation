#' Traceplot of the posterior chain of the proportion (i.e., pi) of unspliced (U) counts in each group - single-cell RNA-seq data
#'
#' \code{plot_traceplot} plots the traceplot of the posterior chain 
#' of the proportion (i.e., pi) of unspliced (U) counts in each group.
#' The vertical grey dashed line indicates the burn-in 
#' (the iterations on the left side of the burn-in are discarded in posterior analyses).
#' 
#' @param results a \code{list} of \code{\linkS4class{data.frame}} objects, 
#' computed via \code{\link{DifferentialRegulation}} (single-cell RNA-seq data), 
#' or \code{\link{DifferentialRegulation_bulk}} (bulk RNA-seq data).
#' @param gene_id a character, indicating the gene to plot.
#' @param cluster_id a character, indicating the cell cluster to plot.
#' 
#' @return A \code{gtable} object.
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
plot_traceplot = function(results,
                          gene_id,
                          cluster_id){
  if(! "MCMC_U" %in% names(results)){
    message("'names(results)' must include 'MCMC_U'.")
    message("Set 'traceplot = TRUE' when running 'DifferentialRegulation' function.")
    return(NULL)
  }
  
  if( !is.list(results) ){
    message("'results' must be a 'list' created via 'DifferentialRegulation' function.")
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
  
  sel_cluster = which( names(results$MCMC_U) == cluster_id)
  
  if( length(sel_cluster) == 0 ){
    message("'cluster_id' not found in 'names(results$MCMC_U)'")
    return(NULL)
  }
  
  sel_gene = which( results$MCMC_U[[sel_cluster]]$Gene_id == gene_id)
  
  if( length(sel_gene) == 0 ){
    message("'gene_id' not found in 'results$MCMC_U[[sel_cluster]]$Gene_id'")
    return(NULL)
  }
  
  group_names = names(results$MCMC_U[[sel_cluster]])[1:2]
  n_iter = nrow(results$MCMC_U[[sel_cluster]][[1]])
  
  DF = data.frame(pi_U_A = results$MCMC_U[[sel_cluster]][[1]][,sel_gene], 
                  pi_U_B = results$MCMC_U[[sel_cluster]][[2]][,sel_gene], 
                  MCMC_iterations = 1:n_iter)
  
  # plot vertical line at burn-in
  sel_cluster_convergence = which( results$Convergence_results$Cluster_id == cluster_id)
  burn_in = results$Convergence_results$burn_in[sel_cluster_convergence]
  
  # Plot the estimated average proportions of each groups:
  ggp1 = ggplot() +
    geom_line(data = DF, aes_string(x = "MCMC_iterations", y = "pi_U_A")) +
    theme_bw() + 
    ylim(0, 1) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size=16), 
          axis.title = element_text(size=14, face="bold"), 
          plot.title = element_text(size=16)) +
    ggtitle(paste("Gene:", gene_id, "- Cluster:", cluster_id, " - ", group_names[1])) +
    ylab(expression(pi[U])) + 
    geom_vline(xintercept = burn_in, linetype="dashed", 
               colour = "darkgrey")

  ggp2 = ggplot() +
    geom_line(data = DF, aes_string(x = "MCMC_iterations", y = "pi_U_B")) +
    theme_bw() + 
    ylim(0, 1) +
    theme(axis.text.x = element_text(vjust = 0.5), 
          axis.text = element_text(size=16), 
          axis.title = element_text(size=14, face="bold"), 
          plot.title = element_text(size=16)) +
    ggtitle(paste("Gene:", gene_id, "- Cluster:", cluster_id, " - ", group_names[2])) +
    xlab("MCMC iteration") +
    ylab(expression(pi[U])) + 
    geom_vline(xintercept = burn_in, linetype="dashed", 
               colour = "darkgrey")
  
  grid.arrange(ggp1, ggp2, nrow = 2)
}