#' Create a list containing the equivalence classes objects object
#'
#' \code{load_EC} imports the equivalence classes (computed by alevin-fry), and stores them into a list.
#' 
#' @param path_to_EC_counts a vector of length equals to the number of samples: 
#' each element indicates the path to the equivalence classes counts of the respective sample (i.e., geqc_counts.mtx file).
#' @param path_to_EC a vector of length equals to the number of samples: 
#' each element indicates the path to the equivalence classes of the respective sample (i.e., gene_eqclass.txt.gz file).
#' @param path_to_cell_id a vector of length equals to the number of samples: 
#' each element indicates the path to the cell ids of the respective sample (i.e., quants_mat_rows.txt file).
#' @param path_to_gene_id a vector of length equals to the number of samples: 
#' each element indicates the path to the gene ids of the respective sample (i.e., quants_mat_cols.txt file).
#' @param sample_ids a vector of length equals to the number of samples: 
#' each element indicates the name of the sample.
#' 
#' @return A \code{list} object.
#' 
#' @examples
#' TODO
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{load_USA}}, \code{\link{DifferentialRegulation}}
#' 
#' @export
load_EC = function(path_to_EC_counts,
                   path_to_EC,
                   path_to_cell_id,
                   path_to_gene_id,
                   sample_ids){
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
  
  # TODO: move to parLapply
  
  EC_counts = list_EC_gene_id = list_EC_USA_id = gene_id_list = list()
  for(i in seq_len(n_samples)){
    # load cell id:
    cell_id = fread(path_to_cell_id[[i]], sep = " ", 
                    quote = "", header = FALSE)[[1]]
    # add sample id in front of cell id to distinguish cells across samples.
    cell_id = paste(sample_ids[i], cell_id, sep = ".")
    
    # load EC counts:
    # ECs should have integers only, as.integer requires less memory than as.numeric
    aa = readMM(path_to_EC_counts[[i]])
    
    # check if counts are sparse matrix: if not, turn counts intro Sparce object:
    if(!is(aa, "dgCMatrix")){
      aa = Matrix(data=aa, 
                  sparse = TRUE)
    }
    EC_counts[[i]] = aa
    # as.integer(aa)
    # this may create issues later!
    
    rownames(EC_counts[[i]]) = cell_id
    
    # load gene EC:
    EC_genes = fread(path_to_EC[[i]],
                     sep = " ", quote = "", header = FALSE)
    
    # gene_id names (corresponding to ids in EC_genes)
    gene_id = fread(paste0(sample_ids[i],"/alevin/quants_mat_cols.txt"), sep = " ", 
                    quote = "", header = FALSE)
    
    n_genes = length(gene_id$V1)/3
    gene_id_list[[i]] = gene_id
    
    EC_USA_type = c( rep("S", n_genes),  rep("U", n_genes), rep("A", n_genes) )
    
    # NOW GENE names are 3 * longer!!!
    # transform genes ids from 1...3*n_genes to 1...n_genes.
    transform_ids = as.integer(rep( 0:(n_genes-1), 3))
    
    n_EC  = length(EC_genes$V1)-2
    # load EC info:
    X = EC_genes$V1[3:(n_EC+2)] # vector with all transcript ids
    # split info separated by "\t":
    X = strsplit(X,"\t",fixed=TRUE)
    # turn character ids into numeric:
    X = lapply(X, as.integer)
    
    EC_id = 1 + sapply(X, function(x){
      x[length(x)]
    }) # corresponds to rows in EC 
    
    EC_gene_id = sapply(X, function(x){
      x[-length(x)] # transform gene ids -> from 1...3*n_genes to 1...n_genes.
    }) # gene ids as in "gene_id" file.
    rm(X)
    
    # ASSOCIATE SUA splicing state, taken from EC_USA_type:
    EC_USA_type = lapply(EC_gene_id, function(x){
      EC_USA_type[x + 1] # I think + 1
    })
    
    # transform gene ids -> from 1...3*n_genes to 1...n_genes. 
    EC_gene_id = sapply(EC_gene_id, function(x){
      transform_ids[x + 1] # transform gene ids -> from 1...3*n_genes to 1...n_genes.
    }) # gene ids as in "gene_id" file.
    n_distinct_genes_per_EC = sapply(EC_gene_id, length)
    
    # turn character ids into numeric:
    EC_USA_type = lapply(EC_USA_type, function(x)
      as.integer(as.integer(factor(x, levels = c("S", "U", "A"))) - 1)
    )
    
    # sort EC so they match col counts in EC_counts
    EC_counts[[i]] = EC_counts[[i]][,EC_id]
    list_EC_gene_id[[i]] = EC_gene_id
    list_EC_USA_id[[i]] = EC_USA_type
    
    # filter ECs with 1 distinct gene only:
    sel_EC = n_distinct_genes_per_EC == 1
    pp = round( 100 * sum(unlist(EC_counts[[i]][,!sel_EC]))/ sum(unlist(EC_counts[[i]])),2) 
    message("The percentage of multi-gene mapping reads in sample '", sample_ids[i], "' is: ", pp)
  }
  
  # check gene_id is the same across samples!
  if(n_samples > 1.5){
    cond = sapply(gene_id_list[2:n_samples], function(x) all(gene_id_list[[1]] == x))
    if(!all(cond) ){
      message("gene ids contained in path_to_gene_id are not the same; maybe samples were not aligned and quantifies with the same reference transcriprome?")
      return(NULL)
    }
  }  
  
  list(EC_counts = EC_counts,
       list_EC_gene_id = list_EC_gene_id,
       list_EC_USA_id = list_EC_USA_id,
       gene_id = gene_id_list[[1]]) # keep only the gene id from the 1st 
}