#' Create a list containing the equivalence classes objects for the bulk RNA-seq data
#'
#' \code{load_bulk_EC} imports the bulk equivalence classes (computed by salmon), and stores them into a list.
#' 
#' @param path_to_eq_classes a vector of length equals to the number of samples: 
#' each element indicates the path to the equivalence classes counts of the respective sample
#' (i.e., aux_info/eq_classes.txt.gz file).
#' @param n_cores the number of cores to parallelize the tasks on.
#' Since parallelization is at the sample level (each sample is parallelized on a thread), 
#' we suggest setting n_cores to the number of sample, as set by default if 'n_cores' is not specified.
#' 
#' @return A \code{list} object.
#' 
#' @examples
#' # load internal data to the package:
#' data_dir = system.file("extdata", package = "DifferentialRegulation")
#' 
#' # specify samples ids:
#' sample_ids = paste0("sample", seq_len(6))
#' 
#' # Equivalence classes:
#' equiv_classes_files = file.path(data_dir, "salmon", sample_ids, "aux_info/eq_classes.txt.gz")
#' file.exists(equiv_classes_files)
#' 
#' # load EC:
#' EC_list = load_bulk_EC(path_to_eq_classes = equiv_classes_files,
#'                        n_cores = 2)
#' 
#' @author Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{load_bulk_US}}, \code{\link{DifferentialRegulation_bulk}}
#' 
#' @export
load_bulk_EC = function(path_to_eq_classes = NULL,
                        n_cores = NULL){ # type can be salmon or kallisto
  if( !all(file.exists(path_to_eq_classes)) ){ # if at least 1 file not found
    message("'path_to_eq_classes' files not found")
    return(NULL)
  }
  
  # define the number of parellel cores (if left unspecified):
  if(is.null(n_cores)){
    message("'n_cores' was left 'NULL'.")
    message("Since tasks are paralellized on samples,
            we will set 'n_cores' to the number of samples.")
    n_cores = length(path_to_eq_classes)
    message("'n_cores' set to: ", n_cores)
  }
  
  # initialize parallel cores (if n_cores > 1)
  if( n_cores > 1.5){ # if n_cores > 1, I use parallel computing tools
    cl = makeCluster(n_cores, setup_strategy = "sequential")
  }else{
    cl = NULL
  }
  
  # load the data:
  x = load_salmon_EC(path_to_eq_classes = path_to_eq_classes,
                     n_cores = n_cores, cl = cl)
  
  EC_counts = lapply(x, function(y){
    y[[1]]
  })
  list_EC_tr_id = lapply(x, function(y){
    y[[2]]
  })
  list_EC_US_id = lapply(x, function(y){
    y[[3]]
  })
  tr_id = x[[1]][[4]] # same for all samples
  
  # report % of counts that are multi-mapping across transcripts (or S/U versions of the same transcript):
  for(i in seq_along(EC_counts)){
    n_distinct_transcripts_per_EC = vapply(list_EC_tr_id[[i]], length, FUN.VALUE = integer(1))
    sel_EC = n_distinct_transcripts_per_EC == 1
    
    pp = round( 100 * sum(EC_counts[[i]][!sel_EC])/ sum(EC_counts[[i]]) )
    message("The percentage of multi-mapping reads in sample ", i, " is: ", pp)
  }
  
  list(EC_counts = EC_counts,
       list_EC_gene_id = list_EC_tr_id,
       list_EC_USA_id = list_EC_US_id,
       gene_id = tr_id)
}

load_salmon_EC = function(path_to_eq_classes,
                          n_cores, cl){
  if( !all(file.exists(path_to_eq_classes)) ){ # if at least 1 file not found
    message("'path_to_eq_classes' files not found")
    return(NULL)
  }
  
  N = length(path_to_eq_classes)
  
  # MAKE SURE THAT "_" is not present in the list of transcript ids:
  # otherwise crease a longer list of "_".
  
  if( n_cores > 1.5){ # if n_cores > 1, I use parallel computing tools
    x = parLapply(cl = cl, X = path_to_eq_classes, fun = salmon_read_eq_classes)
  }else{
    x = lapply(X = path_to_eq_classes, FUN = salmon_read_eq_classes)
  }
  
  return(x)
}

salmon_read_eq_classes = function(fn){
  fr = fread(fn, sep = " ", quote = "", header = FALSE)
  ids = fr$V1[3:(as.integer(fr$V1[1])+2)] # vector with all transcript ids
  ecs = fr$V1[(as.integer(fr$V1[1])+3):nrow(fr)]
  ecs.s = strsplit(ecs,"\t",fixed=TRUE)
  
  cnt = as.integer(vapply(ecs.s, last, FUN.VALUE = character(1)))
  # as.integer(sapply(ecs.s, last)) # counts for each equiv class
  trans = lapply(ecs.s, function(u) as.integer(u[2:(length(u)-1)]))
  # trans: from 0 to N-1 (like in DR, NOT like in BANDITS).
  
  # identify U transcripts (ending in "-U"):
  len = nchar(ids)
  last = substr(ids, len-1, len)
  sel_U = last == "-U"
  
  # turn T/F indicator into 0/1:
  # 0 for S; 1 for U
  sel_U_01 = ifelse(sel_U, 1, 0)
  
  # select S transcripts only:
  ids_S = ids[!sel_U]
  
  # remove "-U" from U transcripts:
  ids[sel_U] = substr(ids[sel_U], 1, nchar(ids[sel_U])-2)
  
  # make a map from ids to ids_S:
  transform_ids = match(ids, ids_S)
  
  # turn character ids into numeric:
  EC_US_type = lapply(trans, function(x){
    sel_U_01[x + 1]
  })
  # manually check a few to ensure it is correct!
  
  # transform gene ids -> from 1...3*n_genes to 1...n_genes. 
  EC_tr_id = lapply(trans, function(x){
    transform_ids[x + 1] - 1
  }) # gene ids as in "gene_id" file.
  
  list(counts=cnt, EC_tr_id=EC_tr_id, EC_US_type=EC_US_type, tr_id = ids_S)
}