test_that("load_bulk_US() works faultlessly.", {
  # load internal data to the package:
  data_dir = system.file("extdata", package = "DifferentialRegulation")
  
  # specify samples ids:
  sample_ids = paste0("sample", seq_len(6))
  
  # US estimates:
  quant_files = file.path(data_dir, "salmon", sample_ids, "quant.sf")

  # load US estimated counts:
  SE = load_bulk_US(quant_files,
                    sample_ids)
  
  expect_is(SE, "SummarizedExperiment")
  expect_true(all( SE$sample_id %in% sample_ids))
})