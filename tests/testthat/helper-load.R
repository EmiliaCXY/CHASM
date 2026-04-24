pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
r_files <- list.files(file.path(pkg_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
invisible(lapply(r_files, source))
