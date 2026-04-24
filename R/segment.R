#' Resolve bin metadata for residual matrices used in segmentation
#'
#' Parses genomic bin labels from residual matrix column names, restores the
#' coordinates of synthetic padding bins, and reorders the residual matrix to
#' match the resolved bin metadata.
#'
#' @param residuals Numeric matrix of residual values with bins in columns and
#'   cells in rows. Column names must encode genomic bin coordinates.
#'
#' @return A list containing the reordered residual matrix (`residuals`) and the
#'   corresponding bin metadata table (`bin_info`).
resolve_residual_bin_info <- function(residuals) {
  bin_info <- data.frame(chrom_bin = colnames(residuals), stringsAsFactors = FALSE)
  bin_info$chrom_bin <- gsub("read_depth_sqrt_centered\\.", "", bin_info$chrom_bin)

  bin_info[, c("chrom", "start", "end", "arm")] <- stringr::str_split_fixed(
    bin_info$chrom_bin,
    "-|_",
    4
  )
  bin_info$start <- as.numeric(bin_info$start)
  bin_info$end <- as.numeric(bin_info$end)
  bin_info <- tidyr::drop_na(bin_info)

  for (row_idx in seq_len(nrow(bin_info))) {
    bin_label <- bin_info[row_idx, "chrom_bin"]
    if (grepl("padding", bin_label)) {
      padding_parts <- strsplit(bin_label, "\\.")[[1]]
      chrom_arm_label <- padding_parts[1]
      chrom_name <- gsub("p", "", chrom_arm_label)
      chrom_name <- gsub("q", "", chrom_name)

      chrom_arm <- "None"
      if (grepl("p", chrom_arm_label)) {
        chrom_arm <- "p"
      } else if (grepl("q", chrom_arm_label)) {
        chrom_arm <- "q"
      }

      padding_index <- as.numeric(padding_parts[3])

      if (chrom_arm != "None") {
        last_end <- max(bin_info[(bin_info$chrom == chrom_name) & (bin_info$arm == chrom_arm), ]$end, na.rm = TRUE)
      } else {
        last_end <- max(bin_info[(bin_info$chrom == chrom_name), ]$end, na.rm = TRUE)
      }

      bin_info[row_idx, "chrom"] <- chrom_name
      bin_info[row_idx, "start"] <- last_end + padding_index
    }
  }

  residuals <- residuals[, paste0("read_depth_sqrt_centered.", bin_info$chrom_bin), drop = FALSE]

  list(
    residuals = residuals,
    bin_info = bin_info
  )
}

#' Segment residual copy-number signal with DNAcopy
#'
#' Resolves bin metadata from residual matrix column names, smooths the signal
#' with \pkg{DNAcopy}, and performs segmentation. Common use only requires the
#' residual matrix, while advanced \pkg{DNAcopy} tuning can be supplied through
#' `segmentation_control`.
#'
#' Supported `segmentation_control` fields are `smooth.region`,
#' `outlier.SD.scale`, `trim`, `min.width`, `undo.splits`, and `verbose`.
#'
#' @param residuals Numeric matrix of residual values with bins in columns and
#'   cells in rows. Column names must encode genomic bin coordinates.
#' @param alpha Significance threshold passed to `DNAcopy::segment()`.
#' @param segmentation_control Optional named list of advanced `DNAcopy`
#'   settings.
#'
#' @return The `output` table from the segmented `DNAcopy` object.
segment_residuals <- function(residuals, alpha = 0.005, segmentation_control = list()) {
  if (is.null(colnames(residuals))) {
    stop("segment_residuals: residuals must have column names of bins.")
  }

  segmentation_defaults <- list(
    smooth.region = 10,
    outlier.SD.scale = 2,
    trim = 0.01,
    min.width = 5,
    undo.splits = "prune",
    verbose = 0
  )
  segmentation_control <- utils::modifyList(segmentation_defaults, segmentation_control)

  residual_setup <- resolve_residual_bin_info(residuals)
  residuals <- residual_setup$residuals
  bin_info <- residual_setup$bin_info

  residual_cna <- DNAcopy::CNA(
    t(residuals),
    bin_info$chrom,
    bin_info$start,
    data.type = "logratio",
    sampleid = rownames(residuals)
  )

  smoothed_cna <- DNAcopy::smooth.CNA(
    residual_cna,
    smooth.region = segmentation_control$smooth.region,
    outlier.SD.scale = segmentation_control$outlier.SD.scale,
    trim = segmentation_control$trim
  )

  segmented_cna <- DNAcopy::segment(
    smoothed_cna,
    alpha = alpha,
    min.width = segmentation_control$min.width,
    undo.splits = segmentation_control$undo.splits,
    verbose = segmentation_control$verbose
  )

  segmented_cna$output
}
