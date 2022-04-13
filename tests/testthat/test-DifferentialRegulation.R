test_that("DifferentialRegulation() works faultlessly.", {
  # load internal data to the package:
  data_dir = system.file("extdata", package = "DifferentialRegulation")
  
  # specify 4 samples ids:
  sample_ids = paste0("sample_", seq_len(4))
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
  
  design = data.frame(sample = sample_ids,
                      group = c("A", "A", "B", "B"))
  design
  
  sce$cell_type = "cell"
  
  # sce-based test:
  set.seed(169612)
  results_USA = DifferentialRegulation(sce = sce,
                                       EC_list = NULL,
                                       design =  design,
                                       sample_col_name = "sample",
                                       group_col_name = "group",
                                       sce_cluster_name = "cell_type",
                                       min_cells_per_cluster = 100, 
                                       min_counts_per_gene_per_group = 20,
                                       n_cores = NULL)

  expect_is(results_USA, "data.frame")
})