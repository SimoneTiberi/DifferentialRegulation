test_that("load_bulk_EC() works faultlessly.", {
  # load internal data to the package:
  data_dir = system.file("extdata", package = "DifferentialRegulation")
  
  # specify samples ids:
  sample_ids = paste0("sample", seq_len(6))
  
  # Equivalence classes:
  equiv_classes_files = file.path(data_dir, "salmon", sample_ids, "aux_info/eq_classes.txt.gz")

  # load EC:
  EC_list = load_bulk_EC(path_to_eq_classes = equiv_classes_files,
                         n_cores = 2)
  
  expect_is(EC_list, "list")
  expect_true( length(EC_list) == 4 )
})