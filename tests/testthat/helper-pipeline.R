chasm_pkg_root <- function() normalizePath(file.path(testthat::test_path(), "..", ".."))
regression_fixture_path <- function() file.path(chasm_pkg_root(), "tests", "testthat", "fixtures", "regression_default_pipeline.rds")
# Compact, order-independent form of a call table: one row per cell x bin with the columns that define the calls.
calls_key <- function(x, cols) { x <- as.data.frame(x)[, cols, drop = FALSE]; x[do.call(order, unname(as.list(x))), , drop = FALSE] }

# Reference pipeline used by the regression (snapshot) tests.
# Runs the package's default CHASM pipeline (as in tutorial/run_chasm_example_full_pipeline.R) on a fixed subset of
# the bundled example matrix, with a fixed seed for robust PCA:
#   normalize_depth -> wavelet_transform -> robust_pca -> inv_wavelet_transform -> segment_residuals -> assign_cn_state
#   -> assign_cn_state.chrom (celltype "unknown") -> merge_calls
# Takes ~0.5 min on 1 CPU (300 cells).
#
# It runs in a FRESH R process (callr) with LC_COLLATE=en_US.UTF-8, sourcing R/ from the package root. Reason (found 2026-11-26):
# CHASM's output depends on the collation set-up. normalize_depth() ends with merge(..., by = "bin"), whose row order follows
# the collation; that order sets the cell order of the wavelet matrix, and the randomized robust PCA (rsvd) then gives a
# different background for a different cell order. testthat sets LC_COLLATE=C (and ICU collation), and resetting the locale
# with Sys.setlocale() inside a test does not restore the outside behaviour, so the pipeline runs in a separate process.
# The reference and all revision results used en_US.UTF-8.
run_reference_pipeline <- function(n_cells = 300L, seed = 1234L, alpha = 0.005,
                                   file = "Read_depth_matrix_chr1_50mb_2_5pct_Rep1.rds", backtransform = "V0", compact = TRUE) {
  callr::r(function(pkg_root, n_cells, seed, alpha, file, backtransform) {
    suppressPackageStartupMessages({ library(dplyr) })
    for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)) source(f)
    ext <- file.path(pkg_root, "inst", "extdata")
    rd <- readRDS(file.path(ext, file))
    rc <- readRDS(file.path(ext, sub("\\.rds$", "_chrom.rds", file)))
    ids <- sort(intersect(rownames(rd), rownames(rc)), method = "radix")[seq_len(n_cells)]   # fixed, locale-independent subset
    rd <- rd[ids, , drop = FALSE]; rc <- rc[ids, , drop = FALSE]
    chromosomes <- paste0("chr", c(1:22, "X", "Y"))

    rownames(rd) <- normalize_dnacopy_ids(rownames(rd)); rd$barcode <- rownames(rd)
    bins <- setdiff(colnames(rd), "barcode"); bins <- bins[grepl("^chr", bins)]
    rd <- rd[, c(bins, "barcode")]
    cd <- normalize_depth(rd, bins)
    wt <- wavelet_transform(cd, bins, chromosomes)
    set.seed(seed); rp <- robust_pca(wt$mat.wavelet.transform)
    res <- inv_wavelet_transform(rp$Sparse_Signal, wt$chrom.informed.wavelet)
    bg <- inv_wavelet_transform(rp$Expected_Normal, wt$chrom.informed.wavelet)
    colnames(res) <- colnames(bg) <- wt$bin.info.padded                    # as in the tutorial
    seg <- segment_residuals(res, alpha = alpha)
    cn_bin <- if (identical(backtransform, "V0")) assign_cn_state(cd, bg, seg) else assign_cn_state(cd, bg, seg, backtransform = backtransform)

    pos <- intersect(chromosomes, colnames(rc))
    rc$ID <- normalize_dnacopy_ids(rownames(rc)); rc$celltype <- "unknown"; rownames(rc) <- rc$ID   # IDs as in the wavelet track (else merge_calls finds no match)
    nb <- assign_cn_state.chrom(rc, pos)
    final <- merge_calls(cn_bin, nb)
    list(cn_bin = cn_bin, nb = nb, final = final, collate = Sys.getlocale("LC_COLLATE"),
         cd = as.data.frame(cd)[, c("barcode", "bin", "lib_size", "beta_i_c", "beta_i_hat")])
  }, args = list(pkg_root = chasm_pkg_root(), n_cells = n_cells, seed = seed, alpha = alpha, file = file, backtransform = backtransform),
     libpath = .libPaths(), env = c(callr::rcmd_safe_env(), LC_COLLATE = "en_US.UTF-8", LC_ALL = "", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1"),
     show = FALSE) -> o
  if (!identical(o$collate, "en_US.UTF-8")) testthat::skip(paste("collation en_US.UTF-8 not available (got", o$collate, ")"))
  if (compact) compact_outputs(o) else o
}
# Columns kept in the reference (keeps the fixture small): everything that defines the calls.
compact_outputs <- function(o) {
  keep <- list(cn_bin = c("ID", "variable", "segment_id", "indicator", "cn_state_raw", "cn_state_adj"),   # raw CN = f(observed, expected); counts dropped to keep the file small
               nb = c("ID", "chrom", "p_value_adj", "called_cna", "cn_state_binom"),
               final = c("ID", "variable", "segment_id", "cn_state_final"))
  lapply(setNames(names(keep), names(keep)), function(n) as.data.frame(o[[n]])[, intersect(keep[[n]], names(o[[n]])), drop = FALSE])
}
