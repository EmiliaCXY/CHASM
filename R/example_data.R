.example_data_files <- function() {
  c(
    "1-27Nx_peakfrag_cov_5mb.rds",
    "1-27Nx_peakfrag_cov_5mb_chrom.rds",
    "Read_depth_matrix_chr1_50mb_2_5pct_Rep1.rds",
    "Read_depth_matrix_chr1_50mb_2_5pct_Rep1_chrom.rds",
    "Read_depth_matrix_chr1_full_2_5pct_Rep1.rds",
    "Read_depth_matrix_chr1_full_2_5pct_Rep1_chrom.rds",
    "Read_depth_matrix_chr1_p_2_5pct_Rep1.rds",
    "Read_depth_matrix_chr1_p_2_5pct_Rep1_chrom.rds",
    "read_depth_spikein_cells_2mb_ALL.rds"
  )
}

.locate_example_data_file <- function(file_name) {
  installed_path <- system.file("extdata", file_name, package = "CHASM")

  if (nzchar(installed_path)) {
    return(installed_path)
  }

  candidate_roots <- c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", "..")
  )
  candidate_paths <- file.path(candidate_roots, "inst", "extdata", file_name)
  existing_path <- candidate_paths[file.exists(candidate_paths)]

  if (length(existing_path) > 0) {
    return(normalizePath(existing_path[[1]]))
  }

  ""
}

#' List bundled example data files
#'
#' Returns the example dataset file names shipped with the package.
#'
#' @return A character vector of bundled example file names.
#' @export
list_example_data <- function() {
  .example_data_files()
}

#' Get the path to a bundled example data file
#'
#' Returns the path to one bundled example RDS file. During development, when
#' the package is sourced from the repository instead of installed, the function
#' falls back to the source-tree \code{inst/extdata} path.
#'
#' @param file_name Name of the bundled example file. Use
#'   [list_example_data()] to see available options.
#'
#' @return A length-one character vector containing the path to the example RDS.
#' @export
example_data_path <- function(file_name) {
  if (!is.character(file_name) || length(file_name) != 1 || is.na(file_name)) {
    stop("example_data_path: file_name must be a single character string.")
  }

  available_files <- .example_data_files()
  if (!(file_name %in% available_files)) {
    stop(
      "example_data_path: file_name must be one of: ",
      paste(available_files, collapse = ", ")
    )
  }

  example_path <- .locate_example_data_file(file_name)
  if (nzchar(example_path)) {
    return(example_path)
  }

  stop("example_data_path: could not locate bundled file '", file_name, "'.")
}

#' Get the bundled example read-depth matrix path
#'
#' Convenience wrapper for the chr1_full wavelet example.
#'
#' @return A length-one character vector containing the path to the example RDS.
#' @export
example_read_depth_path <- function() {
  example_data_path("Read_depth_matrix_chr1_full_2_5pct_Rep1.rds")
}
