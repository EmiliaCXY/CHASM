normalize_dnacopy_ids <- function(x) {
  make.names(x, unique = FALSE)
}

#' Quantify the percentage of the genome with altered copy-number calls
#'
#' Summarizes per-cell copy-number output into the total altered genomic length
#' and percent of assayed bins assigned to a non-diploid copy-number state.
#'
#' @param cn_calls Data frame containing per-bin or repeated segment-level
#'   copy-number assignments.
#' @param bin_ids Optional character vector of genomic bin identifiers. If
#'   `NULL`, the function uses the `chrom_bin` column in `cn_calls` when
#'   present, otherwise the `variable` column.
#' @param id_col Name of the cell identifier column in `cn_calls`.
#' @param segment_col Name of the segment identifier column in `cn_calls`.
#' @param cn_state_col Name of the copy-number state column in `cn_calls`. If
#'   `NULL`, the function uses `cn_state_final`, `cn_state_adj`, or
#'   `CN_state.final` in that order.
#' @param diploid_state Copy-number state treated as unaltered.
#'
#' @return A data frame with one row per cell and columns `ID`,
#'   `altered_length`, `genome_size`, and `pct_genome_altered`.
#' @export
summarize_genome_altered <- function(
  cn_calls,
  bin_ids = NULL,
  id_col = "ID",
  segment_col = "segment_id",
  cn_state_col = "cn_state_final",
  diploid_state = 2
) {
  if (!is.data.frame(cn_calls)) {
    stop("summarize_genome_altered: cn_calls must be a data frame.")
  }

  required_cols <- c(id_col, segment_col)
  missing_cols <- setdiff(required_cols, colnames(cn_calls))
  if (length(missing_cols) > 0) {
    stop(
      "summarize_genome_altered: cn_calls is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  cn_state_col <- .resolve_cn_state_col(cn_calls, cn_state_col)

  if (is.null(bin_ids)) {
    if ("chrom_bin" %in% colnames(cn_calls)) {
      bin_ids <- cn_calls$chrom_bin
    } else if ("variable" %in% colnames(cn_calls)) {
      bin_ids <- cn_calls$variable
    } else {
      stop(
        "summarize_genome_altered: supply bin_ids or provide a 'chrom_bin' or 'variable' column in cn_calls."
      )
    }
  }

  genome_bins <- unique(as.character(bin_ids))
  genome_size <- sum(.interval_lengths(genome_bins, "bin_ids"))

  segment_calls <- unique(
    cn_calls[, c(id_col, segment_col, cn_state_col), drop = FALSE]
  )
  colnames(segment_calls) <- c("ID", "segment_id", "cn_state")
  segment_calls$cn_state <- suppressWarnings(as.numeric(segment_calls$cn_state))

  if (anyNA(segment_calls$cn_state)) {
    stop("summarize_genome_altered: copy-number states must be numeric or coercible to numeric.")
  }

  altered_segments <- segment_calls[segment_calls$cn_state != diploid_state, , drop = FALSE]

  all_ids <- unique(as.character(cn_calls[[id_col]]))
  if (nrow(altered_segments) == 0) {
    return(data.frame(
      ID = all_ids,
      altered_length = 0,
      genome_size = genome_size,
      pct_genome_altered = 0,
      stringsAsFactors = FALSE
    ))
  }

  altered_segments$segment_length <- .interval_lengths(altered_segments$segment_id, "segment_col")

  altered_summary <- stats::aggregate(
    altered_segments$segment_length,
    by = list(ID = altered_segments$ID),
    FUN = sum
  )
  colnames(altered_summary)[2] <- "altered_length"

  summary_df <- merge(
    data.frame(ID = all_ids, stringsAsFactors = FALSE),
    altered_summary,
    by = "ID",
    all.x = TRUE
  )
  summary_df$altered_length[is.na(summary_df$altered_length)] <- 0
  summary_df$genome_size <- genome_size
  summary_df$pct_genome_altered <- (summary_df$altered_length / genome_size) * 100

  summary_df
}

.resolve_cn_state_col <- function(cn_calls, cn_state_col) {
  if (is.null(cn_state_col)) {
    candidate_cols <- c("cn_state_final", "cn_state_adj", "CN_state.final")
    cn_state_col <- candidate_cols[candidate_cols %in% colnames(cn_calls)][1]
    if (is.na(cn_state_col)) {
      stop(
        "summarize_genome_altered: cn_state_col was not supplied and no default state column was found."
      )
    }
  } else if (!(cn_state_col %in% colnames(cn_calls))) {
    stop("summarize_genome_altered: cn_state_col is not present in cn_calls.")
  }

  cn_state_col
}

#' Plot a copy-number heatmap with chromosome gap separators
#'
#' Reorders bins across chromosomes, reshapes merged copy-number calls into a
#' cell-by-bin matrix, and plots a heatmap with column gaps between
#' chromosomes.
#'
#' @param cn_calls Data frame of per-bin copy-number calls, such as the output
#'   of [merge_calls()].
#' @param cn_state_col Name of the copy-number state column in `cn_calls`. If
#'   `NULL`, the function uses `cn_state_final`, `cn_state_adj`, or
#'   `CN_state.final` in that order.
#' @param id_col Name of the cell identifier column in `cn_calls`.
#' @param bin_col Name of the genomic bin column in `cn_calls`.
#' @param chromosomes Chromosome ordering to use when arranging bins.
#' @param cluster_rows Whether to cluster cells in the heatmap.
#' @param show_rownames Whether to display row names.
#' @param show_colnames Whether to display column names.
#' @param max_cn Maximum copy-number value to display in the heatmap.
#' @param min_cn Minimum copy-number value to display in the heatmap.
#' @param color Palette passed to [pheatmap::pheatmap()].
#' @param file Optional file path passed to [pheatmap::pheatmap()].
#' @param ... Additional arguments passed to [pheatmap::pheatmap()].
#'
#' @return A list containing the ordered matrix, chromosome gap positions, bin
#'   metadata, and the `pheatmap` return object.
#' @export
plot_cn_heatmap <- function(
  cn_calls,
  cn_state_col = NULL,
  id_col = "ID",
  bin_col = "chrom_bin",
  chromosomes = c(paste0("chr", 1:22), "chrX", "chrY"),
  cluster_rows = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  max_cn = 4,
  min_cn = 0,
  color = grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
  file = NULL,
  ...
) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("plot_cn_heatmap: package 'pheatmap' must be installed to draw the heatmap.")
  }

  heatmap_input <- .prepare_cn_heatmap_input(
    cn_calls = cn_calls,
    cn_state_col = cn_state_col,
    id_col = id_col,
    bin_col = bin_col,
    chromosomes = chromosomes,
    max_cn = max_cn,
    min_cn = min_cn
  )

  heatmap_obj <- pheatmap::pheatmap(
    heatmap_input$matrix,
    cluster_rows = cluster_rows,
    cluster_cols = FALSE,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    gaps_col = heatmap_input$gaps_col,
    color = color,
    file = file,
    ...
  )

  invisible(c(
    heatmap_input,
    list(heatmap = heatmap_obj)
  ))
}

.prepare_cn_heatmap_input <- function(
  cn_calls,
  cn_state_col = NULL,
  id_col = "ID",
  bin_col = "chrom_bin",
  chromosomes = c(paste0("chr", 1:22), "chrX", "chrY"),
  max_cn = 4,
  min_cn = 0
) {
  if (!is.data.frame(cn_calls)) {
    stop(".prepare_cn_heatmap_input: cn_calls must be a data frame.")
  }

  required_cols <- c(id_col, bin_col)
  missing_cols <- setdiff(required_cols, colnames(cn_calls))
  if (length(missing_cols) > 0) {
    stop(
      ".prepare_cn_heatmap_input: cn_calls is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  cn_state_col <- .resolve_cn_state_col(cn_calls, cn_state_col)

  cn_plot <- unique(cn_calls[, c(id_col, bin_col, cn_state_col), drop = FALSE])
  colnames(cn_plot) <- c("ID", "chrom_bin", "cn_state")

  cn_wide <- tidyr::pivot_wider(
    cn_plot,
    names_from = "chrom_bin",
    values_from = "cn_state"
  )
  cn_wide <- as.data.frame(cn_wide, stringsAsFactors = FALSE)
  rownames(cn_wide) <- cn_wide$ID
  cn_wide$ID <- NULL

  bin_df <- data.frame(chrom_bin = colnames(cn_wide), stringsAsFactors = FALSE)
  bin_df[, c("chrom", "start", "end")] <- stringr::str_split_fixed(bin_df$chrom_bin, "-|_", 3)
  bin_df$chrom <- factor(bin_df$chrom, levels = chromosomes)
  bin_df$start <- as.numeric(bin_df$start)
  bin_df$end <- as.numeric(bin_df$end)

  if (anyNA(bin_df$start) || anyNA(bin_df$end)) {
    stop(".prepare_cn_heatmap_input: bin_col must encode intervals like 'chr1_1_200'.")
  }

  bin_df <- bin_df[order(bin_df$chrom, bin_df$start), , drop = FALSE]
  bin_df$order <- seq_len(nrow(bin_df))

  gap_df <- stats::aggregate(
    bin_df$order,
    by = list(chrom = bin_df$chrom),
    FUN = max
  )
  gap_df <- gap_df[!is.na(gap_df$chrom), , drop = FALSE]
  gaps_col <- gap_df$x
  if (length(gaps_col) > 0) {
    gaps_col <- gaps_col[gaps_col < ncol(cn_wide)]
  }

  cn_wide <- cn_wide[, bin_df$chrom_bin, drop = FALSE]
  cn_matrix <- as.matrix(cn_wide)
  storage.mode(cn_matrix) <- "numeric"
  cn_matrix <- pmin(pmax(cn_matrix, min_cn), max_cn)

  list(
    matrix = cn_matrix,
    gaps_col = gaps_col,
    bin_df = bin_df
  )
}

.interval_lengths <- function(interval_ids, arg_name) {
  interval_parts <- stringr::str_split_fixed(as.character(interval_ids), "-|_", 4)
  starts <- suppressWarnings(as.numeric(interval_parts[, 2]))
  ends <- suppressWarnings(as.numeric(interval_parts[, 3]))

  if (anyNA(starts) || anyNA(ends)) {
    stop(
      "summarize_genome_altered: ", arg_name,
      " must encode genomic intervals like 'chr1_1_200' or 'chr1_1_200_p'."
    )
  }

  ends - starts + 1
}
