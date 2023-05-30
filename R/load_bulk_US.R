#' Creates a SummarizedExperiment containing the estimated US counts for the bulk RNA-seq data
#'
#' \code{load_bulk_US} imports the bulk estimated US (Unspliced, and Spliced) 
#' counts (computed by salmon), and stores them into a \code{SummarizedExperiment} object.
#' 
#' @param path_to_quant_files a vector of length equals to the number of samples: 
#' each element indicates the path to the US estimated count matrix of the respective sample (i.e., quant.sf file).
#' @param sample_ids a vector of length equals to the number of samples,
#' indicating the names of the respective samples in `path_to_quant_files`.
#' 
#' @return A \code{SummarizedExperiment} object.
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
#' # load US estimated counts:
#' SE = load_bulk_US(quant_files,
#'                   sample_ids)
#' 
#' @author Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{load_bulk_EC}}, \code{\link{DifferentialRegulation_bulk}}
#' 
#' @export
load_bulk_US = function(path_to_quant_files,
                        sample_ids){ # type can be salmon or kallisto
  txi = tximport(files = path_to_quant_files, type = "salmon", txOut = TRUE)
  
  counts = txi$counts

  eff_len = txi$length
  eff_len = apply(eff_len, 1, median, na.rm = TRUE)
  
  rm(txi);
  
  # Split S and U transcripts
  len = nchar(rownames(counts))
  last = substr(rownames(counts), len-1, len)
  sel_U = last == "-U"
  
  S = counts[ !sel_U, ]
  U = counts[ sel_U, ]
  
  eff_len_S = eff_len[!sel_U]
  eff_len_U = eff_len[sel_U]
  
  rownames(U) = substr(rownames(U), 1, nchar(rownames(U))-2)
  
  rm(counts); rm(len); rm(last); rm(sel_U); rm(eff_len)
  
  # match names (U -> S) in case U transcripts not in S:
  matches = match(rownames(S), rownames(U) )
  
  U = U[matches,]
  eff_len_U = eff_len_U[matches]
  
  # T/F vector to indicate if transcript has U version:
  has_U_tr = !is.na(matches)
  # set to 0 the NA matches: S transcripts without a corresponding U transcript
  U[!has_U_tr,] = 0
  # set the same length: S transcripts without a corresponding U transcript
  eff_len_U[!has_U_tr] = eff_len_S[!has_U_tr]
  # ensure rownames are identical:
  rownames(U) = rownames(S)
  
  # filter: remove transcripts with 0 counts across both S and U
  sel_non_zeros = rowSums(S + U) > 0
  
  S = S[sel_non_zeros,]
  U = U[sel_non_zeros,]
  has_U_tr = has_U_tr[sel_non_zeros]
  eff_len_S = eff_len_S[sel_non_zeros]
  eff_len_U = eff_len_U[sel_non_zeros]
  
  rm(sel_non_zeros)
  
  SE <- SummarizedExperiment(assays=list(spliced = S, 
                                          unspliced = U),
                              rowData = data.frame(transcript_id = rownames(S),
                                                   has_U_tr = has_U_tr,
                                                   eff_len_S = eff_len_S,
                                                   eff_len_U = eff_len_U),
                              colData = data.frame(sample_id = sample_ids))
  
  rm(S); rm(U)
  
  # filter non-detected transcripts
  tot = assays(SE)$spliced + assays(SE)$unspliced
  sel = rowSums(tot) > 0
  SE = SE[sel,]
  
  rm(tot); rm(sel)
  
  SE
}