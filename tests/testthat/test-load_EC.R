test_that("load_EC() works faultlessly.", {
  # load internal data to the package:
  data_dir = system.file("extdata", package = "DifferentialRegulation")
  
  # specify 4 samples ids:
  sample_ids = paste0("sample_", seq_len(4))
  # set directories of each sample input data (obtained via alevin-fry):
  base_dir = file.path(data_dir, "alevin-fry", sample_ids)
  file.exists(base_dir)
  
  # set paths to USA counts, cell id, gene id, EC counts and ECs:
  path_to_counts = file.path(base_dir,"/alevin/quants_mat.mtx")
  path_to_cell_id = file.path(base_dir,"/alevin/quants_mat_rows.txt")
  path_to_gene_id = file.path(base_dir,"/alevin/quants_mat_cols.txt")
  path_to_EC_counts = file.path(base_dir,"/alevin/geqc_counts.mtx")
  path_to_EC = file.path(base_dir,"/alevin/gene_eqclass.txt.gz")
  
  # load EC counts:
  EC_list = load_EC(path_to_EC_counts,
                    path_to_EC,
                    path_to_cell_id,
                    path_to_gene_id,
                    sample_ids)
  
  expect_is(EC_list, "list")
  expect_true( length(EC_list) == 4 )
})