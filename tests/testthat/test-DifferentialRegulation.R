test_that("DifferentialRegulation() works faultlessly.", {
  # load internal data to the package:
  data_dir = system.file("extdata", package = "DifferentialRegulation")
  
  # specify samples ids:
  sample_ids = paste0("organoid", c(1:3, 16:18))
  # set directories of each sample input data (obtained via alevin-fry):
  base_dir = file.path(data_dir, "alevin-fry", sample_ids)
  file.exists(base_dir)
  
  # set paths to USA counts, cell id and gene id:
  path_to_counts = file.path(base_dir,"/alevin/quants_mat.mtx")
  path_to_cell_id = file.path(base_dir,"/alevin/quants_mat_rows.txt")
  path_to_gene_id = file.path(base_dir,"/alevin/quants_mat_cols.txt")
  
  # load USA counts:
  sce = load_USA(path_to_counts,
                 path_to_cell_id,
                 path_to_gene_id,
                 sample_ids)
  
  # define the design of the study:
  design = data.frame(sample = sample_ids,
                      group = c( rep("3 mon", 3), rep("6 mon", 3) ))
  design
  
  # cell types should be assigned to each cell;
  # here we load pre-computed cell types:
  path_to_DF = file.path(data_dir,"DF_cell_types.txt")
  DF_cell_types = read.csv(path_to_DF, sep = "\t", header = TRUE)
  matches = match(colnames(sce), DF_cell_types$cell_id)
  sce$cell_type = DF_cell_types$cell_type[matches]
  
  # sce-based test:
  set.seed(169612)
  results_USA = DifferentialRegulation(sce = sce,
                                       EC_list = NULL,
                                       design =  design,
                                       sample_col_name = "sample",
                                       group_col_name = "group",
                                       sce_cluster_name = "cell_type",
                                       min_cells_per_cluster = 100, 
                                       min_counts_per_gene_per_group = 20)
  
  expect_is(results_USA, "list")
  expect_is(results_USA[[1]], "data.frame")
  expect_is(results_USA[[2]], "data.frame")
  expect_is(results_USA[[3]], "data.frame")
})