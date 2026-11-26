# Regression (snapshot) test: the DEFAULT pipeline must reproduce the calls of the unmodified code base (7af5239).
# New options (e.g. V2 back-transform, second round) must be off by default, so this test must keep passing.
# ~0.5 min on 1 CPU; skipped on CRAN.
test_that("default pipeline reproduces the 7af5239 reference calls", {
  skip_on_cran()
  fx_path <- regression_fixture_path(); skip_if_not(file.exists(fx_path), "reference fixture missing (tests/testthat/fixtures/make_regression_fixture.R)")
  ref <- readRDS(fx_path)$out
  new <- suppressMessages(run_reference_pipeline())
  for (n in names(ref)) {
    expect_identical(sort(names(new[[n]])), sort(names(ref[[n]])), info = paste(n, "columns"))
    expect_equal(nrow(new[[n]]), nrow(ref[[n]]), info = paste(n, "rows"))
  }
  expect_gt(nrow(ref$final), 0)                                                 # merged calls exist (guards against an ID mismatch)
  # numeric content (tolerance for floating point), compared in a fixed row order
  for (n in names(ref)) {
    cols <- intersect(names(ref[[n]]), names(new[[n]]))
    expect_equal(calls_key(new[[n]], cols), calls_key(ref[[n]], cols), tolerance = 1e-8, ignore_attr = TRUE, info = n)
  }
})
