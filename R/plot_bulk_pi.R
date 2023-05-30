#' Plot the estimated proportions of US counts in each group - bulk RNA-seq data
#'
#' \code{plot_bulk_pi} plots the posterior means of the proportions of US counts, in each group.
#' If 'CI' is TRUE, a profile Wald type confidence interval will also be added;
#' the level of the confidence interval is specified by 'CI_level'.
#' 
#' @param results a \code{list} of 2 \code{\linkS4class{data.frame}} objects, computed via \code{\link{DifferentialRegulation}}.
#' @param transcript_id a character, indicating the transcript to plot.
#' @param CI a logical ('TRUE' by default), indicating whether to plot a 
#' profile Wald type confidence interval around the estimated proportions.
#' @param CI_level a numeric between 0 and 1, indicating the level of the confidence interval.
#' @return A \code{ggplot} object.
#' 
#' @examples
#' # see the example of DifferentialRegulation_bulk function:
#' help(DifferentialRegulation_bulk)
#' 
#' @author Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{DifferentialRegulation}}
#' 
#' @export
plot_bulk_pi = function(results,
                        transcript_id,
                        CI = TRUE,
                        CI_level = 0.95){
  if( (!is.list(results)) | (length(results) != 2) ){
    message("'results' must be a 'list' of 2 objects, as created via 'DifferentialRegulation_bulk' function.")
    return(NULL)
  }
  if(!("Differential_results" %in% names(results))){
    message("'names(results)' must include 'Differential_results'.")
    message("'results' must be a 'list' of 2 objects, as created via 'DifferentialRegulation_bulk' function.")
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
  
  if(!is.character(transcript_id)){
    transcript_id = as.character(transcript_id)
  }
  if( length(transcript_id) != 1){
    message("'transcript_id' contains ", length(transcript_id), " values: one transcript only can be specified.")
    return(NULL)
  }
  
  #select DF and cols for pi and SD
  DF = results$Differential_results
  sel_pi = 5:8
  sel_sd = 9:12
  tr_names = factor(c("S", "U"), levels = c("S", "U") )
  
  group_names = substring(colnames( DF )[c(5,7)],6)
  
  sel = which(DF$Transcript_id == transcript_id)
  
  if( length(sel) == 0 ){
    message("'transcript_id' not found in 'results$Differential_results$Transcript_id'")
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
    ggtitle(paste("Transcript:", transcript_id)) +
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