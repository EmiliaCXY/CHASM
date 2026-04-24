options(scipen = 999)

find_package_root <- function(path) {
  current <- normalizePath(path, winslash = "/", mustWork = FALSE)

  if (file.exists(current) && !dir.exists(current)) {
    current <- dirname(current)
  }

  repeat {
    description_path <- file.path(current, "DESCRIPTION")
    if (file.exists(description_path)) {
      description_lines <- readLines(description_path, warn = FALSE)
      if (any(grepl("^Package:", description_lines))) {
        return(current)
      }
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      return(NULL)
    }
    current <- parent
  }
}

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  arg_match <- grep(file_arg, cmd_args, value = TRUE)
  if (length(arg_match) > 0) {
    return(sub(file_arg, "", arg_match[1]))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(sys.frames()[[1]]$ofile)
  }

  NULL
}

load_chasm <- function() {
  package_root <- find_package_root(getwd())

  if (is.null(package_root)) {
    script_path <- get_script_path()
    if (!is.null(script_path)) {
      package_root <- find_package_root(script_path)
    }
  }

  if (is.null(package_root)) {
    stop(
      "Could not locate the CHASM package root. ",
      "Run this script from the package directory or install CHASM first."
    )
  }

  if (requireNamespace("CHASM", quietly = TRUE)) {
    suppressPackageStartupMessages(library(CHASM))
    return(invisible(TRUE))
  }

  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(package_root, export_all = FALSE, quiet = TRUE)
    return(invisible(TRUE))
  }

  stop(
    "CHASM is not installed and pkgload is unavailable. ",
    "Install the package or install pkgload to load it from the package directory."
  )
}

load_chasm()

suppressPackageStartupMessages({
  library(dplyr)
})

# Update these parameters for the manual dataset you want to evaluate.
spikein.chrom <- "chr1"
spikein.size <- "full"
spikein.cnv <- "-1"
spikein.pct <- "1pct"
spikein.rep <- "Rep1"
bin.size <- "2mb"

input_dir <- "/Users/emiliac/Dropbox/Rotations/Zhang_lab/scATAC/OtherData/Lin_et_al_breast_epi_GSE272504_spike_in_V2"

dir_name <- paste0(spikein.chrom, "_", spikein.size, "_", spikein.cnv, "_", bin.size)
sample_name <- paste0(
  spikein.chrom, "_", spikein.size, "_", spikein.cnv, "_", spikein.pct, "_", spikein.rep
)

read_depth_chrom_path <- file.path(
  input_dir,
  paste0(spikein.chrom, "_", spikein.size),
  "CountMatrix",
  paste0("Read_depth_matrix_", dir_name),
  paste0("Read_depth_matrix_", sample_name, "_chrom.csv")
)

message("Reading chromosome-level read-depth matrix: ", read_depth_chrom_path)
read.depth.chrom <- read.csv(read_depth_chrom_path, row.names = 1, check.names = FALSE)
positions <- colnames(read.depth.chrom)

read.depth.chrom$ID <- rownames(read.depth.chrom)
read.depth.chrom$celltype <- "unknown"
rownames(read.depth.chrom) <- read.depth.chrom$ID

message("Running chromosome-level negative-binomial CN calling")
cn_state.nb <- assign_cn_state.chrom(read.depth.chrom, positions)
cn_state.nb <- cn_state.nb %>% filter(chrom != "chrY")

print(
  cn_state.nb %>%
    dplyr::count(chrom, calledCNA, name = "n_calls")
)

results_path <- file.path("tutorial", paste0(sample_name, "_chrom_nb_calls.csv"))
write.csv(cn_state.nb, results_path, row.names = FALSE)
message("Wrote chromosome-level CN calls to: ", results_path)
