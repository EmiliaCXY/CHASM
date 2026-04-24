recenter_svd_normalized_depth <- function(svd.normalized.read.depth, chrom_depth_per_cell) {
  if (!all(c("barcode", "bin", "beta_i_hat") %in% colnames(chrom_depth_per_cell))) {
    stop("recenter_svd_normalized_depth: chrom_depth_per_cell must contain 'barcode', 'bin', 'beta_i_hat'.")
  }

  chrom_depth_per_cell.beta_hat <- tidyr::pivot_wider(
    chrom_depth_per_cell[, c("barcode", "bin", "beta_i_hat")],
    names_from = "bin",
    values_from = "beta_i_hat"
  )
  rownames(chrom_depth_per_cell.beta_hat) <- chrom_depth_per_cell.beta_hat$barcode
  chrom_depth_per_cell.beta_hat <- as.matrix(chrom_depth_per_cell.beta_hat[, -1, drop = FALSE])

  bins <- data.frame(bin = colnames(svd.normalized.read.depth), stringsAsFactors = FALSE)
  colnames(svd.normalized.read.depth) <- bins$bin
  bins <- bins[bins$bin %in% chrom_depth_per_cell$bin, , drop = FALSE]
  bins_vec <- bins$bin

  svd.normalized.read.depth <- svd.normalized.read.depth[, bins_vec, drop = FALSE]
  chrom_depth_per_cell.beta_hat <- chrom_depth_per_cell.beta_hat[, bins_vec, drop = FALSE]

  svd.normalized.read.depth + chrom_depth_per_cell.beta_hat
}

construct_segment_to_bin_dictionary <- function(output, bins) {
  required_cols <- c("chrom", "loc.start", "loc.end", "ID")
  missing_cols <- setdiff(required_cols, colnames(output))
  if (length(missing_cols) > 0) {
    stop(
      "construct_segment_to_bin_dictionary: output missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  output_df <- output[, c("chrom", "loc.start", "loc.end")] %>% dplyr::distinct()
  colnames(output_df) <- c("chrom", "start", "end")
  output_gr <- GenomicRanges::makeGRangesFromDataFrame(output_df)
  output_df$queryHits <- rownames(output_df)
  colnames(output_df) <- c("chrom", "start_output", "end_output", "queryHits")

  bin_df <- as.data.frame(bins)
  colnames(bin_df) <- c("chrom_bin")

  bin_df[, c("chrom", "start", "end", "arm")] <- stringr::str_split_fixed(bin_df$chrom_bin, "-|_", 4)
  bin_df <- bin_df[, c("chrom", "start", "end")]
  bin_df <- bin_df[bin_df$chrom %in% output_df$chrom, , drop = FALSE]
  rownames(bin_df) <- seq_len(nrow(bin_df))
  bin_gr <- GenomicRanges::makeGRangesFromDataFrame(bin_df)
  bin_df$subjectHits <- rownames(bin_df)
  colnames(bin_df) <- c("chrom", "start_bin", "end_bin", "subjectHits")

  overlap <- GenomicRanges::findOverlaps(output_gr, bin_gr, minoverlap = 2)
  overlap_df <- as.data.frame(overlap)

  overlap_df <- merge(overlap_df, output_df, by = "queryHits")
  overlap_df <- merge(overlap_df, bin_df, by = c("subjectHits", "chrom"))
  overlap_df <- overlap_df[, c("chrom", "start_output", "end_output", "start_bin", "end_bin")]

  output.plot <- merge(
    overlap_df,
    output,
    by.x = c("chrom", "start_output", "end_output"),
    by.y = c("chrom", "loc.start", "loc.end")
  )
  output.plot$chrom_bin <- paste0(output.plot$chrom, "_", output.plot$start_bin, "_", output.plot$end_bin)

  segment_output_bin_dict <- output.plot[, c("chrom", "start_output", "end_output", "start_bin", "end_bin", "ID", "chrom_bin")]
  segment_output_bin_dict$segment_output_ID <- paste0(
    segment_output_bin_dict$chrom, "_",
    segment_output_bin_dict$start_output, "_",
    segment_output_bin_dict$end_output
  )
  segment_output_bin_dict$ID <- gsub("\\.", "-", segment_output_bin_dict$ID)

  segment_output_bin_dict
}

assign_cn_state <- function(chrom_depth_per_cell, svd.normalized.read.depth, segment.output) {
  bins <- data.frame(bin = colnames(svd.normalized.read.depth), stringsAsFactors = FALSE)
  bins$bin <- gsub("read_depth_sqrt_centered\\.", "", bins$bin)
  colnames(svd.normalized.read.depth) <- bins$bin
  bins <- bins[bins$bin %in% chrom_depth_per_cell$bin, , drop = FALSE]
  bins_vec <- bins$bin

  svd.normalized.read.depth <- recenter_svd_normalized_depth(
    svd.normalized.read.depth,
    chrom_depth_per_cell
  )

  observed.read.depth <- chrom_depth_per_cell[, c("barcode", "bin", "beta_i_c")]
  observed.read.depth <- tidyr::pivot_wider(
    observed.read.depth,
    names_from = "bin",
    values_from = "beta_i_c"
  )
  observed.read.depth <- as.data.frame(observed.read.depth)
  rownames(observed.read.depth) <- observed.read.depth$barcode
  observed.read.depth <- observed.read.depth[, -1, drop = FALSE]
  observed.read.depth <- observed.read.depth[rownames(svd.normalized.read.depth), bins_vec, drop = FALSE]

  svd.normalized.read.depth.recenter.sq <- svd.normalized.read.depth^2
  svd.normalized.read.depth.recenter.sq <- as.data.frame(svd.normalized.read.depth.recenter.sq)
  svd.normalized.read.depth.recenter.sq$barcode <- rownames(svd.normalized.read.depth.recenter.sq)

  lib_size <- chrom_depth_per_cell %>%
    dplyr::select(barcode, lib_size) %>%
    dplyr::distinct()
  svd.normalized.read.depth.recenter.sq <- merge(
    svd.normalized.read.depth.recenter.sq,
    lib_size,
    by = "barcode"
  )

  ncol_bins <- length(colnames(svd.normalized.read.depth.recenter.sq)) - 1
  svd.normalized.read.depth.recenter.sq.count <- lapply(
    colnames(svd.normalized.read.depth.recenter.sq)[2:ncol_bins],
    function(col) {
      svd.normalized.read.depth.recenter.sq[[col]] * svd.normalized.read.depth.recenter.sq$lib_size
    }
  )
  svd.normalized.read.depth.recenter.sq.count <- as.data.frame(do.call(cbind, svd.normalized.read.depth.recenter.sq.count))
  rownames(svd.normalized.read.depth.recenter.sq.count) <- svd.normalized.read.depth.recenter.sq$barcode
  colnames(svd.normalized.read.depth.recenter.sq.count) <- colnames(svd.normalized.read.depth.recenter.sq)[2:ncol_bins]
  svd.normalized.read.depth.recenter.sq.count$barcode <- rownames(svd.normalized.read.depth.recenter.sq.count)

  observed.read.depth.long <- as.data.frame(observed.read.depth)
  observed.read.depth.long$barcode <- rownames(observed.read.depth.long)
  observed.read.depth.long <- tidyr::pivot_longer(
    observed.read.depth.long,
    cols = -barcode,
    names_to = "bin",
    values_to = "observed.read.depth"
  )

  observed.read.depth <- observed.read.depth^2
  observed.read.depth$barcode <- rownames(observed.read.depth)
  observed.read.depth <- merge(observed.read.depth, lib_size, by = "barcode")
  observed.read.depth.count <- lapply(bins_vec, function(col) {
    observed.read.depth[[col]] * observed.read.depth$lib_size
  })
  observed.read.depth.count <- as.data.frame(do.call(cbind, observed.read.depth.count))
  rownames(observed.read.depth.count) <- observed.read.depth$barcode
  colnames(observed.read.depth.count) <- bins_vec
  observed.read.depth.count$barcode <- rownames(observed.read.depth.count)

  observed.read.depth.count.melt <- tidyr::pivot_longer(
    observed.read.depth.count,
    cols = -barcode,
    names_to = "variable",
    values_to = "observed_count"
  )
  svd.normalized.read.depth.recenter.sq.count.melt <- tidyr::pivot_longer(
    svd.normalized.read.depth.recenter.sq.count,
    cols = -barcode,
    names_to = "variable",
    values_to = "expected_count"
  )

  comp_df <- merge(
    observed.read.depth.count.melt,
    svd.normalized.read.depth.recenter.sq.count.melt,
    by = c("barcode", "variable")
  )

  comp_df$CN <- round((comp_df$observed_count / comp_df$expected_count) * 2)
  comp_df[, c("chr", "start", "end", "arm")] <- stringr::str_split_fixed(comp_df$variable, "-|_", 4)
  comp_df$start <- as.numeric(comp_df$start)
  comp_df$chrom_bin <- comp_df$variable
  comp_df$chrom_bin <- gsub("_p", "", comp_df$chrom_bin)
  comp_df$chrom_bin <- gsub("_q", "", comp_df$chrom_bin)
  comp_df$chrom_bin <- gsub("_cen", "", comp_df$chrom_bin)

  comp_df$barcode <- gsub("\\#", "-", comp_df$barcode)
  comp_df$barcode <- gsub("\\+", "-", comp_df$barcode)
  comp_df$chrom_bin <- gsub("\\-", "\\_", comp_df$chrom_bin)

  segment_output_bin_dict <- construct_segment_to_bin_dictionary(segment.output, bins_vec)

  chrom_depth_with_segment <- merge(
    segment_output_bin_dict[, c("ID", "chrom_bin", "segment_output_ID")],
    comp_df,
    by.x = c("ID", "chrom_bin"),
    by.y = c("barcode", "chrom_bin")
  )

  chrom_depth_with_segment$chrom_name <- stringr::str_split_fixed(
    chrom_depth_with_segment$chrom_bin,
    "-|_",
    3
  )[, 1]

  chrom_num_seg <- chrom_depth_with_segment %>%
    dplyr::group_by(ID, chrom_name) %>%
    dplyr::summarise(n_segments = dplyr::n_distinct(segment_output_ID), .groups = "drop")

  chrom_num_seg.no_seg <- chrom_num_seg %>% dplyr::filter(n_segments == 1)
  chrom_num_seg.no_seg$identifier <- paste0(chrom_num_seg.no_seg$ID, "_", chrom_num_seg.no_seg$chrom_name)

  chrom_num_seg.has_seg <- chrom_num_seg %>% dplyr::filter(n_segments > 1)
  chrom_num_seg.has_seg$identifier <- paste0(chrom_num_seg.has_seg$ID, "_", chrom_num_seg.has_seg$chrom_name)

  chrom_depth_with_segment.no_seg <- chrom_depth_with_segment %>%
    dplyr::mutate(identifier = paste0(ID, "_", chrom_name)) %>%
    dplyr::filter(identifier %in% chrom_num_seg.no_seg$identifier)

  chrom_depth_with_segment.has_seg <- chrom_depth_with_segment %>%
    dplyr::mutate(identifier = paste0(ID, "_", chrom_name)) %>%
    dplyr::filter(identifier %in% chrom_num_seg.has_seg$identifier)

  cn_bin_no_seg <- chrom_depth_with_segment.no_seg %>%
    dplyr::group_by(ID, segment_output_ID) %>%
    dplyr::mutate(
      total_in_group = sum(expected_count),
      weight = expected_count / total_in_group,
      CN_state_raw = sum(CN * weight),
      CN_state_adj = round(CN_state_raw, 0),
      p_value = if (dplyr::n() >= 3) {
        stats::wilcox.test(expected_count, observed_count, paired = TRUE)$p.value
      } else {
        NA_real_
      }
    )

  cn_bin_no_seg$p_value.adj <- stats::p.adjust(cn_bin_no_seg$p_value, method = "BH")
  cn_bin_no_seg$CN_state_adj <- ifelse(
    cn_bin_no_seg$p_value.adj >= 0.05 & cn_bin_no_seg$CN_state_adj != 2,
    2,
    cn_bin_no_seg$CN_state_adj
  )

  cn_bin_has_seg <- chrom_depth_with_segment.has_seg %>%
    dplyr::group_by(ID, segment_output_ID) %>%
    dplyr::mutate(
      total_in_group = sum(expected_count),
      weight = expected_count / total_in_group,
      CN_state_raw = sum(CN * weight),
      CN_state_adj = round(CN_state_raw, 0)
    )
  cn_bin_has_seg$p_value <- NA_real_
  cn_bin_has_seg$p_value.adj <- NA_real_

  rbind(cn_bin_no_seg, cn_bin_has_seg)
}
