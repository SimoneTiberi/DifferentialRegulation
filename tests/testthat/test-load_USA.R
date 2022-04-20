test_that("load_USA() works faultlessly.", {
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
  
  expect_is(sce, "SingleCellExperiment")
  expect_true(all( sce$sample_id %in% sample_ids))
})