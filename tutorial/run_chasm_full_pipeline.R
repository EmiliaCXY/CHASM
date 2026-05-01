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

  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(package_root, export_all = FALSE, quiet = TRUE)
    return(invisible(TRUE))
  }

  if (requireNamespace("CHASM", quietly = TRUE)) {
    suppressPackageStartupMessages(library(CHASM))
    return(invisible(TRUE))
  }

  stop(
    "CHASM is not installed and pkgload is unavailable. ",
    "Install the package or install pkgload to load it from the package directory."
  )
}

setwd('/Users/emiliac/Documents/Rotations/Zhang_lab/Projects/scATACcnv/package')
load_chasm()


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

chromosomes <- paste0("chr", c(1:22, "X", "Y"))

# Update these parameters for the manual dataset you want to evaluate.
spikein.chrom <- "chr1"
spikein.size <- "50mb"
spikein.cnv <- "2"
spikein.pct <- "5pct"
spikein.rep <- "Rep1"
bin.size <- "2mb"

dir_name_1 <- paste0(spikein.chrom, "_", spikein.size)
dir_name_2 <- paste0(spikein.chrom, "_", spikein.size, "_", spikein.cnv, "_", bin.size)
sample_name <- paste0(spikein.chrom, "_", spikein.size, "_", spikein.cnv, "_", spikein.pct, "_", spikein.rep)

message(
  "Running spike-in with chrom: ", spikein.chrom,
  " size: ", spikein.size,
  " cnv: ", spikein.cnv,
  " pct: ", spikein.pct,
  " rep: ", spikein.rep,
  " bin size: ", bin.size
)

simulate_data_dir <- "/Users/emiliac/Dropbox/Rotations/Zhang_lab/scATAC/OtherData/Lin_et_al_breast_epi_GSE272504_spike_in_V2"
output_root <- file.path("tutorial")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
spikein_cells_path <- file.path(
  simulate_data_dir,
  paste0("read_depth_spikein_cells_", bin.size, "_ALL.csv")
)

sanitize_cell_ids <- function(x) {
  getFromNamespace("normalize_dnacopy_ids", "CHASM")(x)
}

read_depth_path <- file.path(
  simulate_data_dir,
  dir_name_1,
  "CountMatrix",
  paste0("Read_depth_matrix_", dir_name_2),
  paste0("Read_depth_matrix_", sample_name, ".csv")
)

message("Reading read-depth matrix: ", read_depth_path)
read.depth <- read.csv(read_depth_path, row.names = 1, check.names = FALSE)
rownames(read.depth) <- sanitize_cell_ids(rownames(read.depth))
read.depth$barcode <- rownames(read.depth)

bins <- setdiff(colnames(read.depth), "barcode")
ki_gl_bins <- bins[grepl("^(KI|GL)", bins)]
chrom_bins <- bins[grepl("^chr", bins)]

read.depth <- read.depth %>% dplyr::select(-dplyr::all_of(ki_gl_bins))
bins <- chrom_bins

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

message("Segmenting residuals...")
colnames(residuals.mat) <- rownames(wavelet.transform$mat.wavelet.transform)
colnames(svd.normalized.read.depth.mat) <- rownames(wavelet.transform$mat.wavelet.transform)
segment.table <- segment_residuals(residuals.mat, alpha = 0.005)

message("Assigning copy number states...")
cn_bin <- assign_cn_state(chrom.depth.per.cell, svd.normalized.read.depth.mat, segment.table)

wavelet_output <- list(
  cn_bin = cn_bin,
  segment.table = segment.table,
  chrom.depth.per.cell = chrom.depth.per.cell,
  residuals.mat = residuals.mat
)

wavelet_output_path <- file.path(output_root, paste0(sample_name, "_output_wavelet.rds"))
# saveRDS(wavelet_output, file = wavelet_output_path)
message("Saved wavelet output to: ", wavelet_output_path)

read_depth_chrom_path <- file.path(
  simulate_data_dir,
  dir_name_1,
  "CountMatrix",
  paste0("Read_depth_matrix_", dir_name_2),
  paste0("Read_depth_matrix_", sample_name, "_chrom.csv")
)

message("Reading chromosome-level read-depth matrix: ", read_depth_chrom_path)
read.depth.chrom <- read.csv(read_depth_chrom_path, row.names = 1, check.names = FALSE)
positions <- colnames(read.depth.chrom)
read.depth.chrom$ID <- sanitize_cell_ids(rownames(read.depth.chrom))
read.depth.chrom$celltype <- "unknown"
rownames(read.depth.chrom) <- read.depth.chrom$ID

cn_state.nb <- assign_cn_state.chrom(read.depth.chrom, positions)

message("Combining segment-level and chromosome-level results...")
# merge_calls <- function(cn_wavelet, cn_nb) {
#   cn_nb$ID <- sub("\\+", "-", cn_nb$ID)
#   cn_merged <- merge(
#     cn_wavelet,
#     cn_nb,
#     by.x = c("ID", "chrom_name"),
#     by.y = c("ID", "chrom"),
#     suffixes = c("_wl", "_nb")
#   )
#   
#   cn_merged$cn_state_final <- ifelse(
#     cn_merged$called_cna == "YES" & !is.na(cn_merged$p_value_adj_wl),
#     cn_merged$cn_state_binom,
#     ifelse(
#       cn_merged$called_cna == "NO" & !is.na(cn_merged$p_value_adj_wl) & cn_merged$p_value_adj_wl < 0.05,
#       2,
#       cn_merged$cn_state_adj
#     )
#   )
#   cn_merged$cn_state_final <- round(cn_merged$cn_state_final, 0)
#   
#   cn_merged
# }
cn_df.comb <- merge_calls(cn_bin, cn_state.nb)

combined_output_path <- file.path(output_root, paste0(sample_name, "_output_cn_df.comb.rds"))
# saveRDS(cn_df.comb, file = combined_output_path)
message("Saved combined CN calls to: ", combined_output_path)

cn_df.comb$ID_clean <- sanitize_cell_ids(cn_df.comb$ID)

if (file.exists(spikein_cells_path)) {
  message("Reading spike-in cell annotations: ", spikein_cells_path)
  cells.spikein <- read.csv(spikein_cells_path, row.names = 1, check.names = FALSE)
  spikein_ids <- sanitize_cell_ids(rownames(cells.spikein))

  cn_bin.spikein <- cn_df.comb %>%
    filter(ID_clean %in% spikein_ids, chrom_name == spikein.chrom)

  cn_bin.nonspikein <- cn_df.comb %>%
    filter(!ID_clean %in% spikein_ids)

  cn_bin.spikein.false <- cn_df.comb %>%
    filter(!ID_clean %in% spikein_ids, cn_state_final != 2) %>%
    group_by(ID_clean, chrom_name) %>%
    summarise(n = n(), .groups = "drop")

  print(
    cn_bin.spikein %>%
      dplyr::count(cn_state_adj, name = "n_spikein_segments")
  )
  print(
    cn_bin.nonspikein %>%
      dplyr::count(cn_state_adj, name = "n_nonspikein_segments")
  )
  print(head(cn_bin.spikein.false, 10))

  cn_bin.wide <- cn_df.comb %>%
    mutate(ID_plot = ID_clean) %>%
    ungroup() %>%
    select(ID_plot, chrom_bin, cn_state_final) %>%
    pivot_wider(names_from = chrom_bin, values_from = cn_state_final)

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

    heatmap_path <- file.path(output_root, paste0(sample_name, "_heatmap.png"))
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
