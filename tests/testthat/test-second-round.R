# chasm_second_round(): runs end to end on the bundled example (300 cells), returns the documented structure, and is reproducible:
#   (1) the same call twice gives identical final calls;
#   (2) supplying round 1 (from the first run) instead of computing it gives identical final calls (each segmentation is seeded internally).
test_that("chasm_second_round runs and is reproducible", {
  skip_on_cran()
  o <- callr::r(function(pkg_root) {
    suppressPackageStartupMessages({ library(dplyr) })
    for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)) source(f)
    ext <- file.path(pkg_root, "inst", "extdata")
    rd <- readRDS(file.path(ext, "Read_depth_matrix_chr1_50mb_2_5pct_Rep1.rds")); rc <- readRDS(file.path(ext, "Read_depth_matrix_chr1_50mb_2_5pct_Rep1_chrom.rds"))
    ids <- sort(intersect(rownames(rd), rownames(rc)), method = "radix")[1:300]; rd <- rd[ids, ]; rc <- rc[ids, ]
    chromosomes <- paste0("chr", c(1:22, "X", "Y"))
    rownames(rd) <- normalize_dnacopy_ids(rownames(rd)); rd$barcode <- rownames(rd); bins <- setdiff(colnames(rd), "barcode")
    cd <- normalize_depth(rd[, c(bins, "barcode")], bins); wt <- wavelet_transform(cd, bins, chromosomes); rp <- robust_pca(wt$mat.wavelet.transform, seed = 1234)
    pos <- intersect(chromosomes, colnames(rc)); rc$ID <- normalize_dnacopy_ids(rownames(rc)); rc$celltype <- "unknown"; rownames(rc) <- rc$ID
    nb <- assign_cn_state.chrom(rc, pos)
    a <- suppressMessages(chasm_second_round(cd, wt, rp, nb, chromosomes = chromosomes))
    b <- suppressMessages(chasm_second_round(cd, wt, rp, nb, chromosomes = chromosomes))
    c2 <- suppressMessages(chasm_second_round(cd, wt, rp, nb, chromosomes = chromosomes, round1_calls = a$round1))
    list(a = a, b = b, c = c2)
  }, args = list(pkg_root = chasm_pkg_root()), libpath = .libPaths(),
     env = c(callr::rcmd_safe_env(), LC_COLLATE = "en_US.UTF-8", LC_ALL = "", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1"))
  key <- function(x) { x <- as.data.frame(x)[, c("ID", "variable", "segment_id", "cn_state_final")]; x[do.call(order, unname(as.list(x))), ] }
  expect_named(o$a, c("final", "round1", "strict", "candidates", "info"))
  expect_true(all(c("ID", "variable", "segment_id", "cn_state_final") %in% names(o$a$final)))
  expect_true(all(c("cand", "seeds", "threshold", "supported") %in% names(o$a$candidates)))
  expect_gt(nrow(o$a$final), 0)
  expect_equal(key(o$b$final), key(o$a$final), ignore_attr = TRUE)          # same call twice
  expect_equal(key(o$c$final), key(o$a$final), ignore_attr = TRUE)          # round 1 supplied vs computed
  expect_equal(key(o$c$strict), key(o$a$strict), ignore_attr = TRUE)
})
