#' Traceplot of the posterior chain of the proportion (i.e., pi) of unspliced (U) counts in each group - bulk RNA-seq data
#'
#' \code{plot_bulk_traceplot} plots the traceplot of the posterior chain 
#' of the proportion (i.e., pi) of unspliced (U) counts in each group.
#' The vertical grey dashed line indicates the burn-in 
#' (the iterations on the left side of the burn-in are discarded in posterior analyses).
#' 
#' @param results a \code{list} of \code{\linkS4class{data.frame}} objects, 
#' computed via \code{\link{DifferentialRegulation}} (single-cell RNA-seq data), 
#' or \code{\link{DifferentialRegulation_bulk}} (bulk RNA-seq data).
#' @param transcript_id a character, indicating the transcript to plot.
#' 
#' @return A \code{gtable} object.
#' 
#' @examples
#' # see the example of DifferentialRegulation_bulk function:
#' help(DifferentialRegulation_bulk)
#' 
#' @author Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{DifferentialRegulation_bulk}}
#' 
#' @export
plot_bulk_traceplot = function(results,
                               transcript_id){
  if(! "MCMC_U" %in% names(results)){
    message("'names(results)' must include 'MCMC_U'.")
    message("Set 'traceplot = TRUE' when running 'DifferentialRegulation_bulk' function.")
    return(NULL)
  }
  
  if( (!is.list(results)) ){
    message("'results' must be a 'list' created via 'DifferentialRegulation_bulk' function.")
    return(NULL)
  }
  
  if(!is.character(transcript_id)){
    transcript_id = as.character(transcript_id)
  }
  if( length(transcript_id) != 1){
    message("'transcript_id' contains ", length(transcript_id), " values: one gene only can be specified.")
    return(NULL)
  }
  
  sel = which(results$MCMC_U$Transcript_id == transcript_id)
  
  if( length(sel) == 0 ){
    message("'transcript_id' not found in 'results$MCMC_U$Transcript_id'")
    return(NULL)
  }
  
  group_names = names(results$MCMC_U)[1:2]
  n_iter = nrow(results$MCMC_U[[1]])
  
  DF = data.frame(pi_U_A = results$MCMC_U[[1]][,sel], 
                  pi_U_B = results$MCMC_U[[2]][,sel], 
                  MCMC_iterations = 1:n_iter)
  
  # plot vertical line at burn-in
  burn_in = results$Convergence_results$burn_in
  
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
    ggtitle(paste("Transcript:", transcript_id, " - ", group_names[1])) +
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
    ggtitle(paste("Transcript:", transcript_id, " - ", group_names[2])) +
    xlab("MCMC iteration") +
    ylab(expression(pi[U])) + 
    geom_vline(xintercept = burn_in, linetype="dashed", 
               colour = "darkgrey")
  
  grid.arrange(ggp1, ggp2, nrow = 2)
}