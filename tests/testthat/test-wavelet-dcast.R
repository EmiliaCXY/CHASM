# wavelet_transform() now pivots long -> wide with data.table::dcast instead of stats::reshape (speed-up).
# Its output must be identical to the former reshape-based implementation, including the cell (column) order,
# because the randomized robust PCA downstream depends on that order.
reshape_wide_reference <- function(chrom_depth_per_cell) {          # the code removed from wavelet_transform() at 7af5239
  long <- chrom_depth_per_cell[, c("barcode", "bin", "read_depth_sqrt_centered")]
  wide <- stats::reshape(long, idvar = "barcode", timevar = "bin", direction = "wide")
  rownames(wide) <- wide$barcode
  as.matrix(subset(wide, select = -barcode))
}

test_that("dcast pivot in wavelet_transform matches the former reshape pivot exactly", {
  ext <- file.path(chasm_pkg_root(), "inst", "extdata")
  rd <- readRDS(file.path(ext, "Read_depth_matrix_chr1_50mb_2_5pct_Rep1.rds"))
  rd <- rd[sort(rownames(rd), method = "radix")[1:120], , drop = FALSE]
  rownames(rd) <- normalize_dnacopy_ids(rownames(rd)); rd$barcode <- rownames(rd)
  bins <- setdiff(colnames(rd), "barcode")
  cd <- suppressMessages(normalize_depth(rd[, c(bins, "barcode")], bins))
  cd <- cd[sample(nrow(cd)), ]                                        # arbitrary long-table order: cell order must follow first appearance

  wt <- suppressMessages(wavelet_transform(cd, bins))
  ref <- reshape_wide_reference(cd)
  # same cells in the same order as reshape (columns of the wavelet matrix)
  expect_identical(colnames(wt$mat.wavelet.transform), rownames(ref))
  # same values: inverse transform of the coefficients reproduces the reshape matrix (orthogonal wavelet basis)
  back <- inv_wavelet_transform(wt$mat.wavelet.transform, wt$chrom.informed.wavelet)
  colnames(back) <- wt$bin.info.padded
  real <- !grepl("padding", wt$bin.info.padded)
  expect_equal(back[rownames(ref), wt$bin.info.padded[real]], ref[, wt$bin.info.padded[real]], tolerance = 1e-12)
})
