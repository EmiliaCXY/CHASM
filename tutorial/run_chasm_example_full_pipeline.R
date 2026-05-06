rm(list = ls())
options(scipen = 999)

setwd('/Users/emiliac/Documents/Rotations/Zhang_lab/Projects/scATACcnv/package/')

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

resolve_package_root <- function() {
  package_root <- find_package_root(getwd())

  if (is.null(package_root)) {
    script_path <- get_script_path()
    if (!is.null(script_path)) {
      package_root <- find_package_root(script_path)
    }
  }

  package_root
}

load_chasm <- function(package_root) {
  if (!is.null(package_root) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(package_root, export_all = FALSE, quiet = TRUE)
    return(invisible(TRUE))
  }

  if (requireNamespace("CHASM", quietly = TRUE)) {
    suppressPackageStartupMessages(library(CHASM))
    return(invisible(TRUE))
  }

  stop(
    "CHASM is not installed and pkgload is unavailable. ",
    "Run this script from the package directory or install CHASM/pkgload first."
  )
}

sanitize_cell_ids <- function(x) {
  getFromNamespace("normalize_dnacopy_ids", "CHASM")(x)
}

package_root <- resolve_package_root()
load_chasm(package_root)

# ==============================================================================
# Read input data
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

chromosomes <- paste0("chr", c(1:22, "X", "Y"))

# Choose one bundled example read-depth matrix from list_example_data().
read_depth_file <- "Read_depth_matrix_chr1_50mb_2_5pct_Rep1.rds"
read_depth_chrom_file <- sub("\\.rds$", "_chrom.rds", read_depth_file)
spikein_cells_file <- "read_depth_spikein_cells_2mb_ALL.rds"

sample_name <- sub("^Read_depth_matrix_", "", sub("\\.rds$", "", read_depth_file))
spikein.chrom <- sub("^Read_depth_matrix_([^_]+)_.*$", "\\1", read_depth_file)
output_root <- if (!is.null(package_root) && dir.exists(file.path(package_root, "tutorial"))) {
  file.path(package_root, "tutorial")
} else {
  getwd()
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

read_depth_path <- example_data_path(read_depth_file)
read_depth_chrom_path <- example_data_path(read_depth_chrom_file)
spikein_cells_path <- example_data_path(spikein_cells_file)

message("Reading bundled read-depth matrix: ", read_depth_path)
read.depth <- readRDS(read_depth_path)

# ==============================================================================
# Prepare input
# ==============================================================================

rownames(read.depth) <- sanitize_cell_ids(rownames(read.depth))
read.depth$barcode <- rownames(read.depth)

bins <- setdiff(colnames(read.depth), "barcode")
chrom_bins <- bins[grepl("^chr", bins)]
ki_gl_bins <- bins[grepl("^(KI|GL)", bins)]
bins <- setdiff(chrom_bins, ki_gl_bins)
read.depth <- read.depth[, c(bins, "barcode")]

# ==============================================================================
# Normalize and compute cell-specific control
# ==============================================================================

message("Normalizing read depth with wavelet transform and robust PCA...")
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

# colnames(residuals.mat) <- rownames(wavelet.transform$mat.wavelet.transform)
# colnames(svd.normalized.read.depth.mat) <- rownames(wavelet.transform$mat.wavelet.transform)
colnames(residuals.mat) <- wavelet.transform$bin.info.padded
colnames(svd.normalized.read.depth.mat) <- wavelet.transform$bin.info.padded

# ==============================================================================
# Segmentation
# ==============================================================================

message("Segmenting residuals...")
segment.table <- segment_residuals(residuals.mat, alpha = 0.001)


# ==============================================================================
# Assign copy number states
# ==============================================================================

message("Assigning segment-level copy-number states...")
cn_bin <- assign_cn_state(chrom.depth.per.cell, svd.normalized.read.depth.mat, segment.table)

wavelet_output <- list(
  cn_bin = cn_bin,
  segment.table = segment.table,
  chrom.depth.per.cell = chrom.depth.per.cell,
  residuals.mat = residuals.mat
)
wavelet_output_path <- file.path(output_root, paste0(sample_name, "_output_wavelet.rds"))
saveRDS(wavelet_output, file = wavelet_output_path)
message("Saved wavelet output to: ", wavelet_output_path)

# ==============================================================================
# Read input data for chromosomal level analysis
# ==============================================================================

message("Reading bundled chromosome-level matrix: ", read_depth_chrom_path)
read.depth.chrom <- readRDS(read_depth_chrom_path)

# ==============================================================================
# Prepare input
# ==============================================================================

positions <- colnames(read.depth.chrom)
read.depth.chrom$ID <- sanitize_cell_ids(rownames(read.depth.chrom))
read.depth.chrom$celltype <- "unknown" # here we didn't use cell type information in order to compare with other methods 
rownames(read.depth.chrom) <- read.depth.chrom$ID

# ==============================================================================
# Chromosomal level copy number state assignment
# ==============================================================================
message("Assigning chromosome-level copy-number states...")
cn_state.nb <- assign_cn_state.chrom(read.depth.chrom, positions)

# ==============================================================================
# Combine segment-level and chromosome-level results
# ==============================================================================

message("Combining segment-level and chromosome-level results...")
cn_df.comb <- merge_calls(cn_bin, cn_state.nb)
combined_output_path <- file.path(output_root, paste0(sample_name, "_output_cn_df.comb.rds"))
saveRDS(cn_df.comb, file = combined_output_path)
message("Saved combined CN calls to: ", combined_output_path)

# ==============================================================================
# Post-processing
# ==============================================================================

cn_df.comb$ID_clean <- sanitize_cell_ids(cn_df.comb$ID)
cn_df.comb$chrom_name <- str_split_fixed(cn_df.comb$chrom_bin, "_", 3)[, 1]

# ==============================================================================
# Reading spike-in annotations and preparing for plotting
# ==============================================================================

message("Reading bundled spike-in annotations: ", spikein_cells_path)
cells.spikein <- readRDS(spikein_cells_path)
rownames(cells.spikein) <- sub('_', 'X', rownames(cells.spikein))
rownames(cells.spikein) <- sub('-', 'X', rownames(cells.spikein))
spikein_ids <- sanitize_cell_ids(rownames(cells.spikein))


cn_df.plot <- cn_df.comb %>%
  mutate(
    ID_plot = ID_clean,
    spikein_group = ifelse(ID_clean %in% spikein_ids, 0L, 1L)
  ) %>%
  arrange(spikein_group, ID_plot, chrom_name, chrom_bin)

cn_bin.wide <- cn_df.plot %>%
  select(ID_plot, chrom_bin, cn_state_final) %>%
  distinct() %>%
  pivot_wider(names_from = chrom_bin, values_from = cn_state_final) %>% 
  as.data.frame()
rownames(cn_bin.wide) <- cn_bin.wide$ID_plot
cn_bin.wide$ID_plot <- NULL

# organize bins by chromosome and position
bin_df <- data.frame(chrom_bin = colnames(cn_bin.wide), stringsAsFactors = FALSE)
bin_df[, c("chr", "start", "end")] <- stringr::str_split_fixed(bin_df$chrom_bin, "_", 3)
bin_df$chr <- factor(bin_df$chr, levels = chromosomes)
bin_df$start <- as.numeric(bin_df$start)
bin_df <- bin_df %>% dplyr::arrange(chr, start)
bin_df$order <- seq_len(nrow(bin_df))
bin_df_break <- bin_df %>% dplyr::group_by(chr) %>% dplyr::summarise(max_order = max(order), .groups = "drop")

# re-arrange columns in cn_bin.wide according to bin_df
cn_bin.wide <- cn_bin.wide[, bin_df$chrom_bin, drop = FALSE]
cn_bin.wide <- pmin(as.matrix(cn_bin.wide), 4)

if (requireNamespace("pheatmap", quietly = TRUE)) {
  annotation_row <- data.frame(
    Spike_in = ifelse(rownames(cn_bin.wide) %in% spikein_ids, "Yes", "No"),
    row.names = rownames(cn_bin.wide)
  )

  heatmap_path <- file.path(output_root, paste0(sample_name, "_example_heatmap.png"))
  pheatmap::pheatmap(
    cn_bin.wide,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_row = annotation_row,
    annotation_colors = list(Spike_in = c(Yes = "red", No = "blue")),
    gaps_col = bin_df_break$max_order,
    gaps_row = sum(rownames(cn_bin.wide) %in% spikein_ids),
    file = heatmap_path
  )
  message("Wrote heatmap to: ", heatmap_path)
} else {
  message("Skipping heatmap because pheatmap is not installed.")
}

results_path <- file.path(output_root, paste0(sample_name, "_example_cn_calls.csv"))
write.csv(cn_df.comb, results_path, row.names = FALSE)
message("Wrote CN calls to: ", results_path)
