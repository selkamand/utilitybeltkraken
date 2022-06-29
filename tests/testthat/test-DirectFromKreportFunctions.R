test_that("kraken_report_get_child_taxids works", {
  testkreport_path = system.file(package = "utilitybeltkraken", "testfiles/DRR002014_1.100sequences.minikraken8GB.kreport")

  # Expect function to run without error
  expect_error(
    kraken_report_get_child_taxids(testkreport_path, 91347),
    NA
  )

  expect_equal(
    object = kraken_report_get_child_taxids(testkreport_path, 91347, inclusive = TRUE),
    expected = c(91347, 543, 561, 562)
    )

  expect_equal(
    object = kraken_report_get_child_taxids(testkreport_path, 91347, inclusive = FALSE),
    expected = c(543, 561, 562)
  )

  expect_equal(
    object = kraken_report_get_child_taxids(testkreport_path, 0, inclusive = TRUE),
    expected = c(0)
  )

  expect_error(
    object = kraken_report_get_child_taxids(testkreport_path, 0, inclusive = FALSE)
  )


  expect_equal(
    object = kraken_report_get_child_taxids(testkreport_path, 562, inclusive = TRUE),
    expected  = c(562)
  )

  expect_error(
    object = kraken_report_get_child_taxids(testkreport_path, 562, inclusive = FALSE)
  )

})
