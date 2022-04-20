#' Cell-type
#' 
#' @rdname DF_cell_types
#' @name DF_cell_types
#' @aliases DF_cell_types
#' 
#' @param DF_cell_types a \code{data.frame} containing the matching between cell id (1st column) and cell-types (2nd column).
#' The cell-types were obtained from the original study of Velasco et al. (2019): \url{https://www.nature.com/articles/s41586-019-1289-x}.
#' and are publicly available via `meta_combined.txt` file, at: 
#' \url{https://singlecell.broadinstitute.org/single_cell/study/SCP282/reproducible-brain-organoids#study-download}.
#' 
#' The raw data can be downloaded from \url{https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA531650&o=acc_s%3Aa}.
#' 
#' @examples
#' # load internal data to the package:
#' # data_dir = system.file("extdata", package = "DifferentialRegulation")
#' # specify samples ids:
#' # sample_ids = paste0("organoid", c(1:3, 16:18))
#' # set directories of each sample input data (obtained via alevin-fry):
#' # base_dir = file.path(data_dir, "alevin-fry", sample_ids)
#' # file.exists(base_dir)
#' # load USA sce:
#' # path_to_counts = file.path(base_dir,"/alevin/quants_mat.mtx")
#' # path_to_cell_id = file.path(base_dir,"/alevin/quants_mat_rows.txt")
#' # path_to_gene_id = file.path(base_dir,"/alevin/quants_mat_cols.txt")
#' #
#' # file.exists(path_to_counts)
#' # file.exists(path_to_cell_id)
#' # file.exists(path_to_gene_id)
#' # 
#' # sce = load_USA(path_to_counts,
#' #                path_to_cell_id,
#' #                path_to_gene_id,
#' #                sample_ids)
#' # 
#' # # load pre-computed cell-types, from the original study
#' # md <- read.csv("meta_combined.txt", sep = "\t")
#' # md <- md[grepl("PGP1", md$Batch), ]
#' # md$SEQ <- sapply(md$NAME, FUN = function (x) strsplit(x, split = "_")[[1]][3])
#' # 
#' # # keep organoids 1:3 and 16:18 only
#' # md = md[ md$Organoid %in% c("1", "2", "3", "16", "17", "18"), ]
#' # table(md$Organoid)
#' # md$Organoid = as.numeric(md$Organoid)
#' # 
#' # md$cell_id = paste0("organoid", md$Organoid, ".", md$SEQ)
#' # matches = match(colnames(sce), md$cell_id)
#' # 
#' # DF_cell_types = data.frame(cell_id  = md$cell_id[matches],
#' #                            cell_type = md$CellType[matches],
#' #                            row.names = NULL)
#' # table(DF_cell_types$cell_type)
#' # Cycling      RG
#' #    2399    1094 
#' # 
#' # save(DF_cell_types, file = "DF_cell_types.RData")
#' 
#' # load the data.frame and visualize its top
#' data(DF_cell_types, package = "DifferentialRegulation")
#' head(DF_cell_types)
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
NULL