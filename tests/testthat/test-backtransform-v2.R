# V2 back-transform (assign_cn_state(..., backtransform = "V2")).
# Checked against an independent re-implementation of the formula used in the revision analyses
# (revision_claude call_round(): e = lib * ((sqrt(e0 / lib) - beta_med + beta_mean)^2 + var_b), reference cells = lib >= 75th pct;
#  bin CN = round(2 * obs / e); segment CN = sum(bin CN * e) / sum(e)).
test_that("V2 expected counts and copy numbers match the revision formula", {
  skip_on_cran()
  v0 <- run_reference_pipeline(compact = FALSE)
  v2 <- run_reference_pipeline(backtransform = "V2", compact = FALSE)
  # identical segmentation; only expected counts and CN change
  expect_equal(nrow(v2$cn_bin), nrow(v0$cn_bin))
  k <- function(x) paste(x$ID, x$variable, x$segment_id)
  v0b <- as.data.frame(v0$cn_bin); v2b <- as.data.frame(v2$cn_bin)
  v2b <- v2b[match(k(v0b), k(v2b)), ]
  expect_false(anyNA(v2b$ID))
  expect_equal(v2b$observed_count, v0b$observed_count)

  cd <- v0$cd
  libs <- unique(cd[, c("barcode", "lib_size")]); ref <- cd[cd$barcode %in% libs$barcode[libs$lib_size >= quantile(libs$lib_size, 0.75)], ]
  st <- data.frame(bin = sort(unique(ref$bin)))
  st$med <- as.numeric(tapply(ref$beta_i_c, ref$bin, median)[st$bin]); st$mn <- as.numeric(tapply(ref$beta_i_c, ref$bin, mean)[st$bin])
  st$vb <- as.numeric(tapply(ref$beta_i_c, ref$bin, function(v) mean((v - mean(v))^2))[st$bin])
  lib <- libs$lib_size[match(v0b$ID, gsub("[#+]", "-", libs$barcode))]
  i <- match(v0b$variable, st$bin)
  e2 <- lib * ((sqrt(v0b$expected_count / lib) - st$med[i] + st$mn[i])^2 + st$vb[i])
  expect_equal(v2b$expected_count, e2, tolerance = 1e-10)
  bin_cn <- round(2 * v0b$observed_count / e2)
  g <- paste(v0b$ID, v0b$segment_id)
  raw <- c(tapply(bin_cn * e2, g, sum) / tapply(e2, g, sum))
  expect_equal(v2b$cn_state_raw, unname(raw[g]), tolerance = 1e-10)
  expect_equal(v2b$cn_state_adj, round(unname(raw[g])))
})

test_that("V2 lowers the upward bias of raw copy number on unchanged chromosomes", {
  skip_on_cran()
  v0 <- as.data.frame(run_reference_pipeline(compact = FALSE)$cn_bin)
  v2 <- as.data.frame(run_reference_pipeline(backtransform = "V2", compact = FALSE)$cn_bin)
  # the example spike-in only alters chr1, so every other single-segment chromosome should be near CN 2
  s0 <- v0[v0$indicator == "single_segment" & v0$chrom_name != "chr1", ]; s2 <- v2[v2$indicator == "single_segment" & v2$chrom_name != "chr1", ]
  expect_lt(abs(median(s2$cn_state_raw) - 2), abs(median(s0$cn_state_raw) - 2))
})
