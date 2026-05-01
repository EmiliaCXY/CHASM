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

setwd('/Users/emiliac/Documents/Rotations/Zhang_lab/Projects/scATACcnv/package')

package_root <- resolve_package_root()
load_chasm(package_root)

# ==============================================================================
# Read input data
# ==============================================================================
chromosomes <- paste0("chr", c(1:22, "X", "Y"))
sample_name <- "GSM8403676_P17_peakfrag_cov_2mb"
wavelet_file <- paste0(sample_name, ".rds")
chrom_file <- paste0(sample_name, "_chrom.rds")

output_root <- if (!is.null(package_root) && dir.exists(file.path(package_root, "tutorial"))) {
  file.path(package_root, "tutorial")
} else {
  getwd()
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

wavelet_path <- example_data_path(wavelet_file)
chrom_path <- example_data_path(chrom_file)

message("Reading bundled binned read-depth matrix: ", wavelet_path)
read.depth <- readRDS(wavelet_path)

# ==============================================================================
# Prepare input
# ==============================================================================
original_wavelet_ids <- rownames(read.depth)
sanitized_wavelet_ids <- sanitize_cell_ids(original_wavelet_ids)
wavelet_id_map <- data.frame(
  ID = sanitized_wavelet_ids,
  ID_original = original_wavelet_ids,
  stringsAsFactors = FALSE
)
rownames(read.depth) <- sanitized_wavelet_ids
read.depth$barcode <- rownames(read.depth)

bins <- setdiff(colnames(read.depth), "barcode")
chrom_bins <- bins[grepl("^chr", bins)]
ki_gl_bins <- bins[grepl("^(KI|GL)", bins)]
bins <- setdiff(chrom_bins, ki_gl_bins)

# ==============================================================================
# Normalize and compute cell-specific control
# ==============================================================================

message("Normalizing read depth with wavelet transform and robust PCA...")
chrom.depth.per.cell <- normalize_depth(read.depth, bins, chromosomes)
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

# ==============================================================================
# Segmentation
# ==============================================================================

message("Segmenting residuals...")
segment.table <- segment_residuals(residuals.mat, alpha = 0.01)

# ==============================================================================
# Assign copy number states
# ==============================================================================

message("Assigning segment-level copy-number states...")
cn_wavelet <- assign_cn_state(
  chrom.depth.per.cell,
  svd.normalized.read.depth.mat,
  segment.table
)

wavelet_output <- list(
  cn_bin = cn_wavelet,
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

message("Reading bundled chromosome-level matrix: ", chrom_path)
read.depth.chrom <- readRDS(chrom_path)
positions <- colnames(read.depth.chrom)
positions <- positions[grepl("^chr", positions)]

# ==============================================================================
# Prepare input
# ==============================================================================

# Prefix chromosome-level barcodes so they match the wavelet-layer cell IDs.
read.depth.chrom$ID <- sanitize_cell_ids(rownames(read.depth.chrom))
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
cn_df.comb <- merge_calls(cn_wavelet, cn_state.nb)

# ==============================================================================
# Post-processing and save results
# ==============================================================================

cn_df.comb <- merge(cn_df.comb, wavelet_id_map, by = "ID", all.x = TRUE, sort = FALSE)
cn_df.comb$ID <- ifelse(is.na(cn_df.comb$ID_original), cn_df.comb$ID, cn_df.comb$ID_original)
cn_df.comb$ID_original <- NULL
combined_output_path <- file.path(output_root, paste0(sample_name, "_output_cn_df.comb.rds"))
saveRDS(cn_df.comb, file = combined_output_path)
message("Saved combined CN calls to: ", combined_output_path)

message("Summarizing genome altered per cell...")
genome_summary <- summarize_genome_altered(cn_df.comb)
summary_output_path <- file.path(output_root, paste0(sample_name, "_genome_altered_summary.csv"))
write.csv(genome_summary, summary_output_path, row.names = FALSE)
message("Wrote genome altered summary to: ", summary_output_path)
print(utils::head(genome_summary[order(-genome_summary$pct_genome_altered), ], 10))

message("Filtering out cells with >20% genome altered...")
retained_ids <- genome_summary$ID[genome_summary$pct_genome_altered <= 20]
cn_df.filtered <- cn_df.comb[cn_df.comb$ID %in% retained_ids, , drop = FALSE]
message(
  "Retained ", length(unique(cn_df.filtered$ID)), " of ",
  nrow(genome_summary), " cells after CNA burden filtering."
)

if (requireNamespace("pheatmap", quietly = TRUE)) {
  heatmap_path <- file.path(output_root, paste0(sample_name, "_heatmap.png"))
  plot_cn_heatmap(
    cn_df.filtered,
    cn_state_col = 'cn_state_final',
    cluster_rows = TRUE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    file = heatmap_path
  )
  message("Wrote heatmap to: ", heatmap_path)
} else {
  message("Skipping heatmap because pheatmap is not installed.")
}

results_path <- file.path(output_root, paste0(sample_name, "_cn_calls_filtered.csv"))
write.csv(cn_df.filtered, results_path, row.names = FALSE)
message("Wrote filtered CN calls to: ", results_path)
