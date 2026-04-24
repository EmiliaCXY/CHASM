#' Normalize read depth by library size and center by expected bin signal
#'
#' Converts a per-cell read-depth matrix into long format, normalizes each bin
#' by cell library size, and centers the square-root-transformed depth using the
#' median signal among higher-coverage cells.
#'
#' @param df.depth Data frame with one row per cell, a `barcode` column, and one
#'   column per genomic bin.
#' @param bins Character vector of bin column names to include in the
#'   normalization.
#'
#' @return A long-format data frame containing per-cell, per-bin normalized and
#'   centered depth values along with library-size statistics.
normalize_depth <- function(df.depth, bins) {
  df.depth <- as.data.frame(df.depth)

  if (!all(c("barcode", bins) %in% colnames(df.depth))) {
    stop("normalize_depth: df.depth must contain 'barcode' and all bins columns.")
  }

  df.depth$lib_size <- rowSums(df.depth[, bins, drop = FALSE], na.rm = TRUE)

  df.depth_long <- df.depth[, c(bins, "barcode", "lib_size")]
  chrom_depth_per_cell <- tidyr::pivot_longer(
    df.depth_long,
    cols = dplyr::all_of(bins),
    names_to = "bin",
    values_to = "read_depth"
  )

  chrom_depth_per_cell$read_depth_lib_norm <-
    chrom_depth_per_cell$read_depth / chrom_depth_per_cell$lib_size
  chrom_depth_per_cell$beta_i_c <- sqrt(chrom_depth_per_cell$read_depth_lib_norm)

  lib_size_thres2 <- stats::quantile(df.depth$lib_size, probs = 0.75, na.rm = TRUE)

  beta_i_hat_per_chrom <- chrom_depth_per_cell %>%
    dplyr::filter(lib_size >= lib_size_thres2) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(beta_i_hat = stats::median(beta_i_c, na.rm = TRUE), .groups = "drop")

  chrom_depth_per_cell <- merge(
    chrom_depth_per_cell,
    beta_i_hat_per_chrom,
    by = "bin",
    all.x = TRUE
  )

  chrom_depth_per_cell$read_depth_sqrt_centered <-
    chrom_depth_per_cell$beta_i_c - chrom_depth_per_cell$beta_i_hat

  chrom_depth_per_cell
}

#' Construct a chromosome-informed wavelet basis matrix
#'
#' Groups bins by chromosome arm, pads each group to a dyadic length, and places
#' the corresponding Haar wavelet basis on the block diagonal of a single
#' matrix.
#'
#' @param bin_info Data frame describing ordered bins, including chromosome,
#'   arm, and chromosome ordering metadata.
#'
#' @return A list with the updated bin metadata (`bin_info`), per-arm padded bin
#'   counts (`arm_bin_counts`), and the assembled wavelet basis matrix
#'   (`wavelet_matrix`).
construct_chrom_informed_wavelet_matrix <- function(bin_info) {
  message("Constructing chromosome-informed wavelet matrix...")

  bin_info$chrom_arm.adjust <- ifelse(
    bin_info$chromosomes %in% c("chr21", "chr22", "chrY"),
    bin_info$chromosomes,
    bin_info$chrom_arm
  )

  arm_bin_counts <- bin_info %>%
    dplyr::group_by(chrom_arm.adjust, order) %>%
    dplyr::summarise(bin_count = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(order)

  arm_bin_counts$padded_bin_count <- 2 ^ ceiling(log2(arm_bin_counts$bin_count))
  arm_bin_counts <- arm_bin_counts[arm_bin_counts$padded_bin_count > 2, ]

  total_bins <- sum(arm_bin_counts$padded_bin_count)
  wavelet_matrix <- matrix(0, nrow = total_bins, ncol = total_bins)

  row_start <- 1
  col_start <- 1
  for (i in seq_len(nrow(arm_bin_counts))) {
    padded_arm_bin_count <- arm_bin_counts$padded_bin_count[i]

    haar_matrix <- wavethresh::GenW(padded_arm_bin_count, filter.number = 1, family = "DaubExPhase")
    row_end <- row_start + padded_arm_bin_count - 1
    col_end <- col_start + padded_arm_bin_count - 1

    wavelet_matrix[row_start:row_end, col_start:col_end] <- t(haar_matrix)

    row_start <- row_end + 1
    col_start <- col_end + 1
  }

  list(
    bin_info = bin_info,
    arm_bin_counts = arm_bin_counts,
    wavelet_matrix = wavelet_matrix
  )
}

#' Pad centered read-depth data to match the wavelet basis dimensions
#'
#' Expands each chromosome-arm block of the centered read-depth matrix with
#' zero-filled columns so that its width matches the dyadic block size used in
#' the chromosome-informed wavelet matrix.
#'
#' @param centered_depth_matrix Matrix of centered read-depth values with cells
#'   in rows and ordered bins in columns.
#' @param bin_info Data frame describing the ordered bins and their
#'   chromosome-arm assignments.
#' @param arm_bin_counts Data frame containing the original and padded bin
#'   counts for each chromosome arm.
#'
#' @return A padded centered read-depth matrix whose columns align with the
#'   wavelet basis dimensions.
pad_read_depth_matrix_for_wavelet <- function(centered_depth_matrix, bin_info, arm_bin_counts) {
  message("Padding read depth matrix to match wavelet dimensions...")

  padded_depth_matrix <- matrix(0, nrow = nrow(centered_depth_matrix), ncol = 0)
  rownames(padded_depth_matrix) <- rownames(centered_depth_matrix)

  for (i in seq_len(nrow(arm_bin_counts))) {
    chrom_arm <- arm_bin_counts$chrom_arm.adjust[i]
    padded_arm_bin_count <- arm_bin_counts$padded_bin_count[i]

    if (padded_arm_bin_count <= 2) {
      next
    }

    arm_bins <- bin_info[bin_info$chrom_arm.adjust == chrom_arm, "bins"]
    arm_depth <- centered_depth_matrix[, arm_bins, drop = FALSE]

    padding_bin_count <- padded_arm_bin_count - arm_bin_counts$bin_count[i]
    if (padding_bin_count > 0) {
      padding <- matrix(0, nrow = nrow(centered_depth_matrix), ncol = padding_bin_count)
      colnames(padding) <- paste0(chrom_arm, ".padding.", seq_len(padding_bin_count))
      padded_arm_depth <- cbind(arm_depth, padding)
    } else {
      padded_arm_depth <- arm_depth
    }

    padded_depth_matrix <- cbind(padded_depth_matrix, padded_arm_depth)
  }

  padded_depth_matrix
}

#' Apply a chromosome-informed wavelet transformation to centered read depth
#'
#' Reorders bins by chromosome and arm, constructs a chromosome-aware wavelet
#' basis, pads the centered read-depth matrix to the required block sizes, and
#' projects the data into wavelet space.
#'
#' @param chrom_depth_per_cell Long-format data frame containing `barcode`,
#'   `bin`, and `read_depth_sqrt_centered` columns.
#' @param bins Character vector of genomic bin names used to subset and order the
#'   read-depth data.
#' @param chromosomes Character vector giving the chromosome order to apply
#'   before building the wavelet basis.
#'
#' @return A list with the transformed matrix (`mat.wavelet.transform`) and the
#'   chromosome-informed wavelet basis (`chrom.informed.wavelet`).
wavelet_transform <- function(chrom_depth_per_cell, bins, chromosomes) {
  if (!all(c("barcode", "bin", "read_depth_sqrt_centered") %in% colnames(chrom_depth_per_cell))) {
    stop("wavelet_transform: chrom_depth_per_cell must contain 'barcode', 'bin', 'read_depth_sqrt_centered'.")
  }

  message("Normalizing read depth and centering by expected bin factor...")

  centered_depth_long <- chrom_depth_per_cell[, c("barcode", "bin", "read_depth_sqrt_centered")]
  centered_depth_wide <- stats::reshape(
    centered_depth_long,
    idvar = "barcode",
    timevar = "bin",
    direction = "wide"
  )
  rownames(centered_depth_wide) <- centered_depth_wide$barcode
  centered_depth_wide <- subset(centered_depth_wide, select = -barcode)
  centered_depth_matrix <- as.matrix(centered_depth_wide)

  message("Reordering bins by chromosome and position...")

  wide_bin_cols <- colnames(centered_depth_wide)
  bin_info <- data.frame(
    bins = gsub("read_depth_sqrt_centered\\.", "", wide_bin_cols),
    stringsAsFactors = FALSE
  )
  bin_info[, c("chrom", "start", "end", "arm")] <- stringr::str_split_fixed(bin_info$bins, "-|_", 4)
  bin_info$chrom_arm <- paste0(bin_info$chrom, bin_info$arm)
  bin_info$start <- as.numeric(bin_info$start)

  message("Merging bin information with chromosome positions...")

  chromosome_order <- data.frame(chromosomes = chromosomes, order = seq_along(chromosomes))
  bin_info <- merge(chromosome_order, bin_info, by.x = "chromosomes", by.y = "chrom")
  bin_info <- bin_info %>% dplyr::arrange(order, chrom_arm, start)
  bin_info$bins <- paste0("read_depth_sqrt_centered.", bin_info$bins)

  centered_depth_matrix <- centered_depth_matrix[, bin_info$bins, drop = FALSE]

  wavelet_setup <- construct_chrom_informed_wavelet_matrix(bin_info)
  bin_info <- wavelet_setup$bin_info
  arm_bin_counts <- wavelet_setup$arm_bin_counts
  wavelet_matrix <- wavelet_setup$wavelet_matrix

  padded_depth_matrix <- pad_read_depth_matrix_for_wavelet(
    centered_depth_matrix,
    bin_info,
    arm_bin_counts
  )

  padded_depth_transposed <- t(padded_depth_matrix)

  message("Performing wavelet transform on read depth matrix...")

  transformed_matrix <- wavelet_matrix %*% padded_depth_transposed
  rownames(transformed_matrix) <- rownames(padded_depth_transposed)

  list(
    mat.wavelet.transform = transformed_matrix,
    chrom.informed.wavelet = wavelet_matrix
  )
}

#' Decompose wavelet-space data into expected and sparse components
#'
#' Runs robust PCA on the wavelet-transformed matrix to separate the low-rank
#' expected signal from sparse residual structure.
#'
#' @param mat.wavelet.transform Numeric matrix returned by
#'   [wavelet_transform()].
#'
#' @return A list with the low-rank expected component (`Expected_Normal`) and
#'   the sparse residual component (`Sparse_Signal`).
robust_pca <- function(mat.wavelet.transform) {
  if (!requireNamespace("rsvd", quietly = TRUE)) {
    stop("robust_pca: package 'rsvd' must be installed to run robust PCA.", call. = FALSE)
  }

  transformed_matrix <- mat.wavelet.transform
  lambda <- 1 / sqrt(max(nrow(transformed_matrix), ncol(transformed_matrix)))
  rpca_fit <- rsvd::rrpca(transformed_matrix, lambda = lambda)

  expected_normal <- rpca_fit$L
  sparse_signal <- rpca_fit$S

  colnames(expected_normal) <- colnames(mat.wavelet.transform)
  colnames(sparse_signal) <- colnames(mat.wavelet.transform)

  list(
    Expected_Normal = expected_normal,
    Sparse_Signal = sparse_signal
  )
}

#' Project a signal matrix back from wavelet space into bin space
#'
#' Multiplies a wavelet-space signal by the transpose of the wavelet basis and
#' returns the reconstructed signal in the original bin ordering.
#'
#' @param signal_matrix Numeric matrix in wavelet space with signals in columns.
#' @param wavelet_matrix Wavelet basis matrix returned by
#'   [wavelet_transform()].
#'
#' @return A reconstructed matrix in bin space with rows corresponding to the
#'   original signals and columns corresponding to bins.
inv_wavelet_transform <- function(signal_matrix, wavelet_matrix) {
  reconstructed_matrix <- t(wavelet_matrix) %*% signal_matrix
  reconstructed_matrix <- t(reconstructed_matrix)

  rownames(reconstructed_matrix) <- colnames(signal_matrix)
  colnames(reconstructed_matrix) <- rownames(wavelet_matrix)

  reconstructed_matrix
}
