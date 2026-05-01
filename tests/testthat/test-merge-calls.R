test_that("merge_calls uses unsuffixed cn_state_adj from wavelet output", {
  cn_wavelet <- data.frame(
    ID = c("cell-1", "cell-1", "cell-2"),
    chrom = c("chr1", "chr2", "chr1"),
    chrom_name = c("chr1", "chr2", "chr1"),
    chrom_bin = c("chr1_1_100", "chr2_1_100", "chr1_1_100"),
    p_value_adj = c(0.2, 0.001, 0.2),
    cn_state_adj = c(3, 1, 4),
    stringsAsFactors = FALSE
  )

  cn_nb <- data.frame(
    ID = c("cell+1", "cell+1", "cell+2"),
    chrom = c("chr1", "chr2", "chr1"),
    p_value_adj = c(0.5, 0.2, 0.001),
    called_cna = c("NO", "NO", "YES"),
    cn_state_binom = c(NA_real_, NA_real_, 5),
    stringsAsFactors = FALSE
  )

  cn_merged <- merge_calls(cn_wavelet, cn_nb)

  expect_equal(cn_merged$cn_state_final, c(3, 2, 5))
})
