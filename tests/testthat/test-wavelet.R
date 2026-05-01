test_that("normalize_depth centers by beta_i_hat", {
  bins <- c("chr1_1_100_p", "chr1_101_200_p", "chr1_201_300_q", "chr1_301_400_q")
  df.depth <- data.frame(
    barcode = c("cell1", "cell2"),
    chr1_1_100_p = c(2000, 1500),
    chr1_101_200_p = c(1000, 1800),
    chr1_201_300_q = c(1200, 1400),
    chr1_301_400_q = c(1800, 1300),
    stringsAsFactors = FALSE
  )

  chrom_depth <- normalize_depth(df.depth, bins)
  expect_true(all(c("read_depth_sqrt_centered", "beta_i_c", "beta_i_hat") %in% colnames(chrom_depth)))

  center_diff <- chrom_depth$beta_i_c - chrom_depth$beta_i_hat
  expect_lt(max(abs(center_diff - chrom_depth$read_depth_sqrt_centered)), 1e-10)
})

test_that("normalize_depth reorders bins by chromosome and position", {
  bins <- c("chr2_101_200_p", "chr1_201_300_q", "chr1_1_100_p")
  df.depth <- data.frame(
    barcode = "cell1",
    chr2_101_200_p = 30,
    chr1_201_300_q = 20,
    chr1_1_100_p = 10,
    stringsAsFactors = FALSE
  )

  chrom_depth <- normalize_depth(df.depth, bins, chromosomes = c("chr1", "chr2"))

  expect_equal(
    chrom_depth$bin,
    c("chr1_1_100_p", "chr1_201_300_q", "chr2_101_200_p")
  )
})

test_that("wavelet_transform returns consistent dimensions", {
  bins <- c("chr1_1_100_p", "chr1_101_200_p", "chr1_201_300_q", "chr1_301_400_q")
  df.depth <- data.frame(
    barcode = c("cell1", "cell2"),
    chr1_1_100_p = c(2000, 1500),
    chr1_101_200_p = c(1000, 1800),
    chr1_201_300_q = c(1200, 1400),
    chr1_301_400_q = c(1800, 1300),
    stringsAsFactors = FALSE
  )
  chrom_depth <- normalize_depth(df.depth, bins)

  wt <- wavelet_transform(chrom_depth, bins, chromosomes = "chr1")

  expect_true(is.matrix(wt$mat.wavelet.transform))
  expect_true(is.matrix(wt$chrom.informed.wavelet))
  expect_equal(nrow(wt$chrom.informed.wavelet), ncol(wt$chrom.informed.wavelet))
  expect_equal(ncol(wt$mat.wavelet.transform), length(unique(chrom_depth$barcode)))
})

test_that("robust_pca decomposes matrix", {
  skip_if_not_installed("rsvd")

  m <- matrix(c(1, 2, 3, 4, 10, 2), nrow = 3)
  colnames(m) <- c("cell1", "cell2")
  rpca_out <- robust_pca(m)

  expect_true(all(c("Expected_Normal", "Sparse_Signal") %in% names(rpca_out)))
  expect_equal(dim(rpca_out$Expected_Normal), dim(m))
  expect_equal(dim(rpca_out$Sparse_Signal), dim(m))
  expect_lt(max(abs(rpca_out$Expected_Normal + rpca_out$Sparse_Signal - m)), 1e-4)
})

test_that("inv_wavelet_transform restores dimensions", {
  mat.signal <- matrix(1:6, nrow = 3)
  colnames(mat.signal) <- c("cell1", "cell2")

  mat.wavelet <- diag(3)
  rownames(mat.wavelet) <- c("b1", "b2", "b3")

  inv <- inv_wavelet_transform(mat.signal, mat.wavelet)
  expect_equal(dim(inv), c(2, 3))
  expect_lt(max(abs(inv - t(mat.signal))), 1e-10)
})

test_that("construct_segment_to_bin_dictionary maps segments to bins", {
  segment_table <- data.frame(
    chrom = c("chr1", "chr1"),
    loc.start = c(1, 201),
    loc.end = c(200, 300),
    ID = c("cell1", "cell1"),
    stringsAsFactors = FALSE
  )
  bins <- c("chr1_1_100_p", "chr1_101_200_p", "chr1_201_300_q")

  segment_bin_map <- construct_segment_to_bin_dictionary(segment_table, bins)

  expect_equal(nrow(segment_bin_map), 3)
  expect_equal(
    unique(segment_bin_map$segment_id),
    c("chr1_1_200", "chr1_201_300")
  )
  expect_equal(
    segment_bin_map$chrom_bin,
    c("chr1_1_100", "chr1_101_200", "chr1_201_300")
  )
})
