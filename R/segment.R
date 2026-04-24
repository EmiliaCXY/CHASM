segment_residuals <- function(residuals, alpha = 0.005) {
  if (is.null(colnames(residuals))) {
    stop("segment_residuals: residuals must have column names of bins.")
  }

  chrom_bins <- data.frame(chrom_bin = colnames(residuals), stringsAsFactors = FALSE)
  chrom_bins$chrom_bin <- gsub("read_depth_sqrt_centered\\.", "", chrom_bins$chrom_bin)

  chrom_bins[, c("chrom", "start", "end", "arm")] <- stringr::str_split_fixed(
    chrom_bins$chrom_bin,
    "-|_",
    4
  )
  chrom_bins$start <- as.numeric(chrom_bins$start)
  chrom_bins$end <- as.numeric(chrom_bins$end)
  chrom_bins <- tidyr::drop_na(chrom_bins)

  for (row in seq_len(nrow(chrom_bins))) {
    loc <- chrom_bins[row, "chrom_bin"]
    if (grepl("padding", loc)) {
      split_text <- strsplit(loc, "\\.")[[1]]
      chr_arm <- split_text[1]
      chr_name <- gsub("p", "", chr_arm)
      chr_name <- gsub("q", "", chr_name)

      arm <- "None"
      if (grepl("p", chr_arm)) {
        arm <- "p"
      } else if (grepl("q", chr_arm)) {
        arm <- "q"
      }

      pad_num <- as.numeric(split_text[3])

      if (arm != "None") {
        last_loc <- max(chrom_bins[(chrom_bins$chrom == chr_name) & (chrom_bins$arm == arm), ]$end, na.rm = TRUE)
      } else {
        last_loc <- max(chrom_bins[(chrom_bins$chrom == chr_name), ]$end, na.rm = TRUE)
      }

      chrom_bins[row, "chrom"] <- chr_name
      chrom_bins[row, "start"] <- last_loc + pad_num
    }
  }

  residuals <- residuals[, paste0("read_depth_sqrt_centered.", chrom_bins$chrom_bin), drop = FALSE]

  residuals.CNA.object <- DNAcopy::CNA(
    t(residuals),
    chrom_bins$chrom,
    chrom_bins$start,
    data.type = "logratio",
    sampleid = rownames(residuals)
  )

  smoothed.CNA.object <- DNAcopy::smooth.CNA(
    residuals.CNA.object,
    smooth.region = 10,
    outlier.SD.scale = 2,
    trim = 0.01
  )

  segment.smoothed.CNA.object <- DNAcopy::segment(
    smoothed.CNA.object,
    alpha = alpha,
    min.width = 5,
    undo.splits = "prune",
    verbose = 0
  )

  segment.smoothed.CNA.object$output
}
