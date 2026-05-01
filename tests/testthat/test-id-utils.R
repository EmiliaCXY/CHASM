test_that("normalize_dnacopy_ids matches DNAcopy-compatible name mangling", {
  raw_ids <- c("A+B", "A-B", "1-27Nx|cell-1", "if", "cell name")

  expect_equal(
    normalize_dnacopy_ids(raw_ids),
    c("A.B", "A.B", "X1.27Nx.cell.1", "if.", "cell.name")
  )
})
