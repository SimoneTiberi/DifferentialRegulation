#' Creates a SingleCellExperiment containing the estimated USA counts for the single-cell RNA-seq data
#'
#' \code{load_USA} imports the single-cell estimated USA (Unspliced, Spliced and Ambiguous) 
#' counts (computed by alevin-fry), and stores them into a \code{SingleCellExperiment} object.
#' 
#' @param path_to_counts a vector of length equals to the number of samples: 
#' each element indicates the path to the USA estimated count matrix of the respective sample (i.e., quants_mat.mtx file).
#' @param path_to_cell_id a vector of length equals to the number of samples: 
#' each element indicates the path to the cell ids of the respective sample (i.e., quants_mat_rows.txt file).
#' @param path_to_gene_id a vector of length equals to the number of samples: 
#' each element indicates the path to the gene ids of the respective sample (i.e., quants_mat_cols.txt file).
#' @param sample_ids a vector of length equals to the number of samples: 
#' each element indicates the name of the sample.
#' 
#' @return A \code{SingleCellExperiment} object.
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
#' @author Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{load_EC}}, \code{\link{DifferentialRegulation}}
#' 
#' @export
load_USA = function(path_to_counts,
                    path_to_cell_id,
                    path_to_gene_id,
                    sample_ids){
  if( !all(file.exists(path_to_counts)) ){ # if at least 1 file not found
    message("'path_to_counts' files not found")
    return(NULL)
  }
  if( !all(file.exists(path_to_cell_id)) ){ # if at least 1 file not found
    message("'path_to_cell_id' files not found")
    return(NULL)
  }
  if( !all(file.exists(path_to_gene_id)) ){ # if at least 1 file not found
    message("'path_to_gene_id' files not found")
    return(NULL)
  }
  
  n_samples = length(sample_ids)
  
  if( n_samples != length(path_to_counts) ){
    message("sample_ids, path_to_counts, path_to_cell_id and path_to_gene_id must have the same length")
    return(NULL)
  }
  if( n_samples != length(path_to_cell_id) ){
    message("sample_ids, path_to_counts, path_to_cell_id and path_to_gene_id must have the same length")
    return(NULL)
  }
  if( n_samples != length(path_to_gene_id) ){
    message("sample_ids, path_to_counts, path_to_cell_id and path_to_gene_id must have the same length")
    return(NULL)
  }
  
  spliced = list()
  unspliced = list()
  ambiguous = list()
  
  # TODO: move to parLapply
  
  for(i in seq_len(n_samples)){
    counts = readMM(path_to_counts[[i]])
    n_genes = ncol(counts)/3
    
    # first n_genes refer to S
    # second n_genes refer to U
    # third n_genes refer to A
    spliced[[i]]  = counts[, seq.int(1, n_genes) ]
    unspliced[[i]] = counts[, seq.int(n_genes+1, 2*n_genes)]
    ambiguous[[i]] = counts[, seq.int(2*n_genes+1, 3*n_genes)]
    
    # load cell id:
    cell_id = fread(path_to_cell_id[[i]], sep = " ", 
                    quote = "", header = FALSE)[[1]]
    # add sample id in front of cell id to distinguish cells across samples.
    cell_id = paste(sample_ids[i], cell_id, sep = ".")
    
    # load gene_id id:
    gene_id = fread(path_to_gene_id[[i]], sep = " ", 
                    quote = "", header = FALSE)[[1]]
    gene_id = gene_id[ seq.int(1, n_genes) ]
    
    rownames(spliced[[i]]) = cell_id
    colnames(spliced[[i]]) = gene_id
    rownames(unspliced[[i]]) = cell_id
    colnames(unspliced[[i]]) = gene_id
    rownames(ambiguous[[i]]) = cell_id
    colnames(ambiguous[[i]]) = gene_id
  }
  
  # ensure the length of spliced is the same across samples (i.e., same number of genes):
  if(length(unique(vapply(spliced, ncol, FUN.VALUE = integer(1) ))) != 1){
    message("counts in path_to_counts must have the same number of genes")
    return(NULL)
  }
  
  S = t(do.call(rbind, spliced))
  U = t(do.call(rbind, unspliced))
  A = t(do.call(rbind, ambiguous))
  
  # we define "counts" as spliced + 50% of ambiguous counts
  C = S + 0.5 * A
  
  sce <- SingleCellExperiment(assays=list(spliced= S, 
                                          unspliced= U,
                                          ambiguous= A,
                                          counts = C),
                              colData = data.frame(sample_id = factor(rep(sample_ids, vapply(spliced, nrow,  FUN.VALUE = integer(1) ))) ))
  # check that rownames = gene id
  # check that colnames = cell id
  
  sce
}
