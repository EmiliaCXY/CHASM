#' Recenter normalized depth by the expected bin-specific signal
#'
#' Aligns the SVD-normalized depth matrix with the bins present in the
#' normalized long-format depth table and adds back the per-bin expected signal
#' estimate (`beta_i_hat`).
#'
#' @param svd_normalized_depth Numeric matrix of normalized depth values with
#'   cells in rows and bins in columns.
#' @param chrom_depth_per_cell Long-format data frame containing `barcode`,
#'   `bin`, and `beta_i_hat` columns.
#'
#' @return A recentered depth matrix with the same dimensions as the aligned
#'   input matrix.
recenter_svd_normalized_depth <- function(svd_normalized_depth, chrom_depth_per_cell) {
  if (!all(c("barcode", "bin", "beta_i_hat") %in% colnames(chrom_depth_per_cell))) {
    stop("recenter_svd_normalized_depth: chrom_depth_per_cell must contain 'barcode', 'bin', 'beta_i_hat'.")
  }

  beta_hat_matrix <- tidyr::pivot_wider(
    chrom_depth_per_cell[, c("barcode", "bin", "beta_i_hat")],
    names_from = "bin",
    values_from = "beta_i_hat"
  )
  rownames(beta_hat_matrix) <- beta_hat_matrix$barcode
  beta_hat_matrix <- as.matrix(beta_hat_matrix[, -1, drop = FALSE])

  bin_info <- data.frame(bin = colnames(svd_normalized_depth), stringsAsFactors = FALSE)
  colnames(svd_normalized_depth) <- bin_info$bin
  bin_info <- bin_info[bin_info$bin %in% chrom_depth_per_cell$bin, , drop = FALSE]
  bin_names <- bin_info$bin

  svd_normalized_depth <- svd_normalized_depth[, bin_names, drop = FALSE]
  beta_hat_matrix <- beta_hat_matrix[, bin_names, drop = FALSE]

  svd_normalized_depth + beta_hat_matrix
}

#' Map segmentation intervals back to overlapping genomic bins
#'
#' Builds genomic ranges for segment calls and bin coordinates, finds their
#' overlaps, and returns a table linking each segment to the bins it spans.
#'
#' @param segment_table Data frame of segmentation calls containing `chrom`,
#'   `loc.start`, `loc.end`, and `ID` columns.
#' @param bins Character vector of genomic bin labels.
#'
#' @return A data frame mapping each segment call to overlapping bins, including
#'   a derived segment identifier (`segment_id`).
construct_segment_to_bin_dictionary <- function(segment_table, bins) {
  required_cols <- c("chrom", "loc.start", "loc.end", "ID")
  missing_cols <- setdiff(required_cols, colnames(segment_table))
  if (length(missing_cols) > 0) {
    stop(
      "construct_segment_to_bin_dictionary: output missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  segment_ranges <- segment_table[, c("chrom", "loc.start", "loc.end")] %>% dplyr::distinct()
  segment_ranges <- as.data.frame(segment_ranges, stringsAsFactors = FALSE)
  colnames(segment_ranges) <- c("chrom", "start_output", "end_output")
  segment_ranges$start_output <- as.numeric(segment_ranges$start_output)
  segment_ranges$end_output <- as.numeric(segment_ranges$end_output)

  if (anyNA(segment_ranges$start_output) || anyNA(segment_ranges$end_output)) {
    stop("construct_segment_to_bin_dictionary: segment_table contains non-numeric loc.start/loc.end values.")
  }

  segment_granges <- GenomicRanges::GRanges(
    seqnames = segment_ranges$chrom,
    ranges = IRanges::IRanges(
      start = segment_ranges$start_output,
      end = segment_ranges$end_output
    )
  )

  bin_ranges <- as.data.frame(bins, stringsAsFactors = FALSE)
  colnames(bin_ranges) <- c("chrom_bin")

  bin_ranges[, c("chrom", "start", "end", "arm")] <- stringr::str_split_fixed(bin_ranges$chrom_bin, "-|_", 4)
  bin_ranges <- bin_ranges[, c("chrom_bin", "chrom", "start", "end")]
  bin_ranges <- bin_ranges[bin_ranges$chrom %in% segment_ranges$chrom, , drop = FALSE]
  colnames(bin_ranges) <- c("chrom_bin", "chrom", "start_bin", "end_bin")
  bin_ranges$start_bin <- as.numeric(bin_ranges$start_bin)
  bin_ranges$end_bin <- as.numeric(bin_ranges$end_bin)

  if (anyNA(bin_ranges$start_bin) || anyNA(bin_ranges$end_bin)) {
    stop(
      "construct_segment_to_bin_dictionary: bins must encode numeric start/end coordinates, e.g. 'chr1_1_2000000_p'."
    )
  }

  bin_granges <- GenomicRanges::GRanges(
    seqnames = bin_ranges$chrom,
    ranges = IRanges::IRanges(
      start = bin_ranges$start_bin,
      end = bin_ranges$end_bin
    ),
    chrom_bin = bin_ranges$chrom_bin
  )

  segment_bin_overlap <- GenomicRanges::findOverlaps(segment_granges, bin_granges, minoverlap = 2)
  
  segment_ranges_df <- as.data.frame(segment_granges)[S4Vectors::queryHits(segment_bin_overlap), , drop = FALSE]
  colnames(segment_ranges_df) <- c("chrom", "start_output", "end_output", "width_output", 'strand_output')
  bin_granges_df <- as.data.frame(bin_granges)[S4Vectors::subjectHits(segment_bin_overlap), , drop = FALSE]
  colnames(bin_granges_df) <- c("chrom2", "start_bin", "end_bin", "width_bin", 'strand_bin', 'chrom_bin')
  overlap_table <- cbind(segment_ranges_df, bin_granges_df)
  overlap_table <- overlap_table[, c("chrom", "start_output", "end_output", "start_bin", "end_bin")]
  rownames(overlap_table) <- NULL

  segment_bin_table <- merge(
    overlap_table,
    segment_table,
    by.x = c("chrom", "start_output", "end_output"),
    by.y = c("chrom", "loc.start", "loc.end")
  )
  segment_bin_table$chrom_bin <- paste0(segment_bin_table$chrom, "_", segment_bin_table$start_bin, "_", segment_bin_table$end_bin)

  segment_bin_map <- segment_bin_table[, c("chrom", "start_output", "end_output", "start_bin", "end_bin", "ID", "chrom_bin")]
  segment_bin_map$segment_id <- paste0(
    segment_bin_map$chrom, "_",
    segment_bin_map$start_output, "_",
    segment_bin_map$end_output
  )
  segment_bin_map$ID <- gsub("\\.", "-", segment_bin_map$ID)

  segment_bin_map
}

#' Convert normalized depth values into expected or observed read counts
#'
#' Squares the depth matrix, joins per-cell library sizes, and scales each bin
#' by library size to recover count-like values for downstream comparison.
#'
#' @param depth_matrix Numeric matrix of depth values with cells in rows and
#'   bins in columns.
#' @param lib_size Data frame containing `barcode` and `lib_size` columns.
#' @param bin_names Character vector giving the bin columns to retain and order.
#'
#' @return A data frame of per-cell, per-bin counts with a `barcode` column.
construct_count_matrix <- function(depth_matrix, lib_size, bin_names) {
  depth_sq <- depth_matrix^2
  depth_sq <- as.data.frame(depth_sq)
  depth_sq$barcode <- rownames(depth_sq)

  depth_sq <- merge(
    depth_sq,
    lib_size,
    by = "barcode"
  )

  count_matrix <- lapply(
    bin_names,
    function(col) {
      depth_sq[[col]] * depth_sq$lib_size
    }
  )
  count_matrix <- as.data.frame(do.call(cbind, count_matrix))
  rownames(count_matrix) <- depth_sq$barcode
  colnames(count_matrix) <- bin_names
  count_matrix$barcode <- rownames(count_matrix)

  count_matrix
}

#' Construct copy-number calls for single-segment and multi-segment chromosomes
#'
#' Splits per-bin count comparisons into chromosomes with a single segment
#' versus multiple segments, then computes weighted copy-number summaries for
#' each segment. Single-segment chromosomes additionally receive Wilcoxon-based
#' p-values and BH-adjusted calls.
#'
#' @param counts_with_segments Data frame linking per-bin observed and expected
#'   counts to segment identifiers and chromosome labels.
#'
#' @return A data frame of segment-level copy-number calls with raw and adjusted
#'   copy-number states.
construct_cn_calls_by_segment_type <- function(counts_with_segments) {
  segment_counts_by_chrom <- counts_with_segments %>%
    dplyr::group_by(ID, chrom_name) %>%
    dplyr::summarise(n_segments = dplyr::n_distinct(segment_id), .groups = "drop")

  single_segment_chroms <- segment_counts_by_chrom %>% dplyr::filter(n_segments == 1)
  single_segment_chroms$cell_chrom_key <- paste0(single_segment_chroms$ID, "_", single_segment_chroms$chrom_name)

  multi_segment_chroms <- segment_counts_by_chrom %>% dplyr::filter(n_segments > 1)
  multi_segment_chroms$cell_chrom_key <- paste0(multi_segment_chroms$ID, "_", multi_segment_chroms$chrom_name)

  single_segment_counts <- counts_with_segments %>%
    dplyr::mutate(cell_chrom_key = paste0(ID, "_", chrom_name)) %>%
    dplyr::filter(cell_chrom_key %in% single_segment_chroms$cell_chrom_key)

  multi_segment_counts <- counts_with_segments %>%
    dplyr::mutate(cell_chrom_key = paste0(ID, "_", chrom_name)) %>%
    dplyr::filter(cell_chrom_key %in% multi_segment_chroms$cell_chrom_key)

  cn_calls_single_segment <- single_segment_counts %>%
    dplyr::group_by(ID, segment_id) %>%
    dplyr::mutate(
      total_in_group = sum(expected_count),
      weight = expected_count / total_in_group,
      cn_state_raw = sum(bin_cn * weight),
      cn_state_adj = round(cn_state_raw, 0),
      p_value = if (dplyr::n() >= 3) {
        stats::wilcox.test(expected_count, observed_count, paired = TRUE)$p.value
      } else {
        NA_real_
      }
    )

  cn_calls_single_segment$p_value_adj <- stats::p.adjust(cn_calls_single_segment$p_value, method = "BH")
  cn_calls_single_segment$cn_state_adj <- ifelse(
    cn_calls_single_segment$p_value_adj >= 0.05 & cn_calls_single_segment$cn_state_adj != 2,
    2,
    cn_calls_single_segment$cn_state_adj
  )

  cn_calls_multi_segment <- multi_segment_counts %>%
    dplyr::group_by(ID, segment_id) %>%
    dplyr::mutate(
      total_in_group = sum(expected_count),
      weight = expected_count / total_in_group,
      cn_state_raw = sum(bin_cn * weight),
      cn_state_adj = round(cn_state_raw, 0)
    )
  cn_calls_multi_segment$p_value <- NA_real_
  cn_calls_multi_segment$p_value_adj <- NA_real_

  rbind(cn_calls_single_segment, cn_calls_multi_segment)
}

#' Assign segment-level copy-number states from observed and expected depth
#'
#' Reconstructs observed and expected read counts, links them to segmentation
#' intervals, and summarizes each segment into a copy-number state.
#'
#' @param chrom_depth_per_cell Long-format normalized depth table returned by
#'   [normalize_depth()].
#' @param svd_normalized_depth Recentered expected depth matrix in bin space.
#' @param segment_table Segmentation output returned by [segment_residuals()].
#'
#' @return A data frame of segment-level copy-number calls.
assign_cn_state <- function(chrom_depth_per_cell, svd_normalized_depth, segment_table) {
  bin_info <- data.frame(bin = colnames(svd_normalized_depth), stringsAsFactors = FALSE)
  bin_info$bin <- gsub("read_depth_sqrt_centered\\.", "", bin_info$bin)
  colnames(svd_normalized_depth) <- bin_info$bin
  bin_info <- bin_info[bin_info$bin %in% chrom_depth_per_cell$bin, , drop = FALSE]
  bin_names <- bin_info$bin

  svd_normalized_depth <- recenter_svd_normalized_depth(
    svd_normalized_depth,
    chrom_depth_per_cell
  )

  observed_depth_matrix <- chrom_depth_per_cell[, c("barcode", "bin", "beta_i_c")]
  observed_depth_matrix <- tidyr::pivot_wider(
    observed_depth_matrix,
    names_from = "bin",
    values_from = "beta_i_c"
  )
  observed_depth_matrix <- as.data.frame(observed_depth_matrix)
  rownames(observed_depth_matrix) <- observed_depth_matrix$barcode
  observed_depth_matrix <- observed_depth_matrix[, -1, drop = FALSE]
  observed_depth_matrix <- observed_depth_matrix[rownames(svd_normalized_depth), bin_names, drop = FALSE]

  lib_size <- chrom_depth_per_cell %>%
    dplyr::select(barcode, lib_size) %>%
    dplyr::distinct()
  expected_count_matrix <- construct_count_matrix(svd_normalized_depth, lib_size, bin_names)

  observed_count_matrix <- construct_count_matrix(observed_depth_matrix, lib_size, bin_names)

  observed_count_long <- tidyr::pivot_longer(
    observed_count_matrix,
    cols = -barcode,
    names_to = "variable",
    values_to = "observed_count"
  )
  expected_count_long <- tidyr::pivot_longer(
    expected_count_matrix,
    cols = -barcode,
    names_to = "variable",
    values_to = "expected_count"
  )

  count_comparison <- merge(
    observed_count_long,
    expected_count_long,
    by = c("barcode", "variable")
  )

  count_comparison$bin_cn <- round((count_comparison$observed_count / count_comparison$expected_count) * 2)
  count_comparison[, c("chrom", "start", "end", "arm")] <- stringr::str_split_fixed(count_comparison$variable, "-|_", 4)
  count_comparison$start <- as.numeric(count_comparison$start)
  count_comparison$chrom_bin <- count_comparison$variable
  count_comparison$chrom_bin <- gsub("_(p|q|cen)$", "", count_comparison$chrom_bin)

  count_comparison$barcode <- gsub("\\#", "-", count_comparison$barcode)
  count_comparison$barcode <- gsub("\\+", "-", count_comparison$barcode)
  count_comparison$chrom_bin <- gsub("\\-", "\\_", count_comparison$chrom_bin)

  segment_bin_map <- construct_segment_to_bin_dictionary(segment_table, bin_names)

  counts_with_segments <- merge(
    segment_bin_map[, c("ID", "chrom_bin", "segment_id")],
    count_comparison,
    by.x = c("ID", "chrom_bin"),
    by.y = c("barcode", "chrom_bin")
  )

  counts_with_segments$chrom_name <- stringr::str_split_fixed(
    counts_with_segments$chrom_bin,
    "-|_",
    3
  )[, 1]

  construct_cn_calls_by_segment_type(counts_with_segments)
}
