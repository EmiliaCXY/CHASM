# Creates tests/testthat/fixtures/regression_default_pipeline.rds: the reference output of the DEFAULT CHASM pipeline
# (run_reference_pipeline() in helper-pipeline.R) on the UNMODIFIED code base (EmiliaCXY/CHASM main @ 7af5239, 2026-05-06).
# Run once from the package root, before any code change:  Rscript tests/testthat/fixtures/make_regression_fixture.R
# Do NOT regenerate after changing R/ unless the change to the default output is intended and reviewed.
suppressPackageStartupMessages({ library(dplyr) })
for (f in list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)) source(f)
test_path <- function(...) file.path("tests", "testthat", ...)
source(file.path("tests", "testthat", "helper-pipeline.R"))
if (requireNamespace("testthat", quietly = TRUE)) assignInNamespace("test_path", test_path, "testthat")
t0 <- Sys.time(); out <- run_reference_pipeline()
commit <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA_character_)
dirty <- tryCatch(length(system("git status --porcelain R", intern = TRUE)) > 0, error = function(e) NA)
fx <- list(out = out, meta = list(created = format(Sys.time()), git_commit = commit, R_dir_modified = dirty, n_cells = 300L, seed = 1234L, alpha = 0.005,
           R = R.version.string, DNAcopy = as.character(packageVersion("DNAcopy")), wavethresh = as.character(packageVersion("wavethresh")),
           minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)))
saveRDS(fx, regression_fixture_path(), compress = "xz")
str(fx$meta); for (n in names(out)) cat(n, ":", paste(dim(out[[n]]), collapse = " x "), "|", paste(names(out[[n]]), collapse = ", "), "\n")
