test_that("estimate_param computes expected_depth", {
  positions <- c("chr1", "chr2", "chr3")
  read_depth <- data.frame(
    ID = c("cell1", "cell2", "cell3", "cell4"),
    celltype = c("A", "A", "B", "B"),
    chr1 = c(2000, 1800, 2200, 2100),
    chr2 = c(2200, 1900, 2000, 2300),
    chr3 = c(1900, 1700, 2100, 2000),
    stringsAsFactors = FALSE
  )

  chrom_depth_per_cell <- estimate_param(read_depth, positions)
  expect_true(all(c("lib_size", "beta_i_hat", "gamma_ci", "expected_depth") %in% colnames(chrom_depth_per_cell)))
  expect_true(all(chrom_depth_per_cell$expected_depth > 0))
  expect_lt(
    max(abs(chrom_depth_per_cell$expected_depth - chrom_depth_per_cell$lib_size * chrom_depth_per_cell$beta_i_hat)),
    1e-10
  )
})

test_that("estimate_variance_dispersion adds phi", {
  positions <- c("chr1", "chr2", "chr3")
  read_depth <- data.frame(
    ID = c("cell1", "cell2", "cell3", "cell4"),
    celltype = c("A", "A", "B", "B"),
    chr1 = c(2000, 1800, 2200, 2100),
    chr2 = c(2200, 1900, 2000, 2300),
    chr3 = c(1900, 1700, 2100, 2000),
    stringsAsFactors = FALSE
  )

  chrom_depth_per_cell <- estimate_param(read_depth, positions)
  chrom_depth_per_cell2 <- estimate_variance_dispersion(chrom_depth_per_cell, coverage_bin_size = 0.01)
  expect_true("phi" %in% colnames(chrom_depth_per_cell2))
  expect_true(any(is.finite(chrom_depth_per_cell2$phi)))
})

test_that("cn_test_nb returns expected columns", {
  positions <- c("chr1", "chr2", "chr3")
  read_depth <- data.frame(
    ID = c("cell1", "cell2", "cell3", "cell4"),
    celltype = c("A", "A", "B", "B"),
    chr1 = c(2000, 1800, 2200, 2100),
    chr2 = c(2200, 1900, 2000, 2300),
    chr3 = c(1900, 1700, 2100, 2000),
    stringsAsFactors = FALSE
  )

  chrom_depth_per_cell <- estimate_param(read_depth, positions)
  chrom_depth_per_cell2 <- estimate_variance_dispersion(chrom_depth_per_cell, coverage_bin_size = 0.01)
  cn_state_df <- cn_test_nb(chrom_depth_per_cell2)
  expect_true(all(c("ID", "chrom", "p_val_adjust", "calledCNA", "CN_state.binom") %in% colnames(cn_state_df)))
  expect_true(all(cn_state_df$p_val_adjust >= 0 & cn_state_df$p_val_adjust <= 1))
  expect_true(all(cn_state_df$calledCNA %in% c("YES", "NO")))
})

test_that("assign_cn_state.chrom returns expected columns", {
  positions <- c("chr1", "chr2", "chr3")
  read_depth <- data.frame(
    ID = c("cell1", "cell2", "cell3", "cell4"),
    celltype = c("A", "A", "B", "B"),
    chr1 = c(2000, 1800, 2200, 2100),
    chr2 = c(2200, 1900, 2000, 2300),
    chr3 = c(1900, 1700, 2100, 2000),
    stringsAsFactors = FALSE
  )

  cn_state_df2 <- assign_cn_state.chrom(read_depth, positions)
  expect_true(all(c("ID", "chrom", "p_val_adjust", "calledCNA", "CN_state.binom") %in% colnames(cn_state_df2)))
})
