rm(list = ls())
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

setwd('/Users/emiliac/Documents/Rotations/Zhang_lab/Projects/scATACcnv/package/')
load_chasm()

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

chromosomes <- paste0("chr", c(1:22, "X", "Y"))

# Update these parameters for the manual dataset you want to evaluate.
spikein.chrom <- "chr1"
spikein.size <- "full"
spikein.cnv <- "2"
spikein.pct <- "10pct"
spikein.rep <- "Rep1"
bin.size <- "2mb"

input_dir <- "/Users/emiliac/Dropbox/Rotations/Zhang_lab/scATAC/OtherData/Lin_et_al_breast_epi_GSE272504_spike_in_V2"

dir_name <- paste0(spikein.chrom, "_", spikein.size, "_", spikein.cnv, "_", bin.size)
sample_name <- paste0(
  spikein.chrom, "_", spikein.size, "_", spikein.cnv, "_", spikein.pct, "_", spikein.rep
)

read_depth_path <- file.path(
  input_dir,
  paste0(spikein.chrom, "_", spikein.size),
  "CountMatrix",
  paste0("Read_depth_matrix_", dir_name),
  paste0("Read_depth_matrix_", sample_name, ".csv")
)

spikein_cells_path <- file.path(
  input_dir,
  paste0("read_depth_spikein_cells_", bin.size, "_ALL.csv")
)

sanitize_cell_ids <- function(x) {
  getFromNamespace("normalize_dnacopy_ids", "CHASM")(x)
}

message("Reading read-depth matrix: ", read_depth_path)
read.depth <- read.csv(read_depth_path, row.names = 1, check.names = FALSE)
rownames(read.depth) <- sanitize_cell_ids(rownames(read.depth))
read.depth$barcode <- rownames(read.depth)

bins <- setdiff(colnames(read.depth), "barcode")
chrom_bins <- bins[grepl("^chr", bins)]
ki_gl_bins <- bins[grepl("^(KI|GL)", bins)]
bins <- setdiff(chrom_bins, ki_gl_bins)
read.depth <- read.depth[, c(bins, "barcode")]

message("Running normalization and wavelet pipeline")
chrom.depth.per.cell <- normalize_depth(read.depth, bins)
wavelet.transform <- wavelet_transform(chrom.depth.per.cell, bins, chromosomes)
robust.pca <- robust_pca(wavelet.transform$mat.wavelet.transform)

residuals.mat <- inv_wavelet_transform(
  robust.pca$Sparse_Signal,
  wavelet.transform$chrom.informed.wavelet
)
svd.normalized.read.depth.mat <- inv_wavelet_transform(
  robust.pca$Expected_Normal,
  wavelet.transform$chrom.informed.wavelet
)

colnames(residuals.mat) <- rownames(wavelet.transform$mat.wavelet.transform)
colnames(svd.normalized.read.depth.mat) <- rownames(wavelet.transform$mat.wavelet.transform)

message("Segmenting residuals and assigning copy-number states")
segment.output <- segment_residuals(residuals.mat, alpha = 0.005)

cn_bin <- assign_cn_state(
  chrom.depth.per.cell,
  svd.normalized.read.depth.mat,
  segment.output
)

cn_bin$ID_clean <- sanitize_cell_ids(cn_bin$ID)

if (file.exists(spikein_cells_path)) {
  message("Reading spike-in cell annotations: ", spikein_cells_path)
  cells.spikein <- read.csv(spikein_cells_path, row.names = 1, check.names = FALSE)
  spikein_ids <- sanitize_cell_ids(rownames(cells.spikein))

  cn_bin.spikein <- cn_bin %>%
    filter(ID_clean %in% spikein_ids, chrom_name == spikein.chrom)

  cn_bin.nonspikein <- cn_bin %>%
    filter(!ID_clean %in% spikein_ids)

  cn_bin.spikein.false <- cn_bin %>%
    filter(!ID_clean %in% spikein_ids, CN_state_adj != 2) %>%
    group_by(ID_clean, chrom_name) %>%
    summarise(n = n(), .groups = "drop")

  print(
    cn_bin.spikein %>%
      dplyr::count(CN_state_adj, name = "n_spikein_segments")
  )
  print(
    cn_bin.nonspikein %>%
      dplyr::count(CN_state_adj, name = "n_nonspikein_segments")
  )
  print(head(cn_bin.spikein.false, 10))

  cn_bin.wide <- cn_bin %>%
    mutate(ID_plot = ID_clean) %>%
    ungroup() %>%
    select(ID_plot, chrom_bin, CN_state_adj) %>%
    pivot_wider(names_from = chrom_bin, values_from = CN_state_adj)

  cn_bin.wide <- as.data.frame(cn_bin.wide)
  rownames(cn_bin.wide) <- cn_bin.wide$ID_plot
  cn_bin.wide$ID_plot <- NULL

  bin_df <- data.frame(chrom_bin = colnames(cn_bin.wide), stringsAsFactors = FALSE)
  bin_df[, c("chr", "start", "end")] <- str_split_fixed(bin_df$chrom_bin, "_", 3)
  bin_df$chr <- factor(bin_df$chr, levels = chromosomes)
  bin_df$start <- as.numeric(bin_df$start)
  bin_df <- bin_df %>% arrange(chr, start)
  bin_df$order <- seq_len(nrow(bin_df))
  bin_df_break <- bin_df %>% group_by(chr) %>% summarise(max_order = max(order), .groups = "drop")

  cn_bin.wide <- cn_bin.wide[, bin_df$chrom_bin, drop = FALSE]
  cn_bin.wide <- cn_bin.wide[c(unique(cn_bin.spikein$ID_clean), unique(cn_bin.nonspikein$ID_clean)), , drop = FALSE]
  cn_bin.wide <- pmin(as.matrix(cn_bin.wide), 4)

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    annotation_row <- data.frame(
      Spike_in = ifelse(rownames(cn_bin.wide) %in% spikein_ids, "Yes", "No"),
      row.names = rownames(cn_bin.wide)
    )

    heatmap_path <- file.path("tutorial", paste0(sample_name, "_heatmap.png"))
    pheatmap::pheatmap(
      cn_bin.wide,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = FALSE,
      show_colnames = FALSE,
      annotation_row = annotation_row,
      annotation_colors = list(Spike_in = c(Yes = "red", No = "blue")),
      gaps_col = bin_df_break$max_order,
      gaps_row = length(unique(cn_bin.spikein$ID_clean)),
      file = heatmap_path
    )
    message("Wrote heatmap to: ", heatmap_path)
  } else {
    message("Skipping heatmap because pheatmap is not installed.")
  }
} else {
  message("Spike-in annotation file not found; skipping spike-in validation and heatmap.")
}

results_path <- file.path("tutorial", paste0(sample_name, "_cn_calls.csv"))
write.csv(cn_bin, results_path, row.names = FALSE)
message("Wrote CN calls to: ", results_path)
