#' Estimate per-cell and per-region parameters for chromosome-level CN testing
#'
#' Computes per-cell library sizes, filters cells by library-size thresholds,
#' estimates region-specific expected depth factors, and derives normalized
#' copy-number ratios and residuals for chromosome-level analysis.
#'
#' @param read_depth Data frame containing chromosome-level read depths along
#'   with `ID` and `celltype` columns.
#' @param positions Character vector of chromosome-level count columns to use.
#' @param min_lib_size Minimum library size required for a cell to be retained.
#' @param max_lib_size Maximum library size allowed for a cell to be retained.
#'
#' @return A long-format data frame with per-cell, per-chromosome expected depth
#'   and copy-number parameters.
estimate_param <- function(read_depth, positions, min_lib_size = 5000, max_lib_size = 60000) {
  message("Step1: Estimating cell specific factor: library size")

  read_depth <- as.data.frame(read_depth)
  if (!all(c("celltype", "ID") %in% colnames(read_depth))) {
    stop("estimate_param: read_depth must contain 'celltype' and 'ID' columns.")
  }
  if (!all(positions %in% colnames(read_depth))) {
    stop("estimate_param: positions not found in read_depth columns.")
  }

  read_depth$lib_size <- rowSums(read_depth[, positions, drop = FALSE], na.rm = TRUE)
  read_depth <- read_depth[read_depth$lib_size >= min_lib_size, , drop = FALSE]
  read_depth <- read_depth[read_depth$lib_size <= max_lib_size, , drop = FALSE]

  message("Step2: Estimating region specific factor")

  chrom_depth_per_cell <- tidyr::pivot_longer(
    read_depth,
    cols = dplyr::all_of(positions),
    names_to = "chrom",
    values_to = "read_depth"
  )

  chrom_depth_per_cell$beta_i_c <- chrom_depth_per_cell$read_depth / chrom_depth_per_cell$lib_size

  beta_i_hat_per_chrom <- chrom_depth_per_cell %>%
    dplyr::group_by(chrom, celltype) %>%
    dplyr::summarise(beta_i_hat = stats::median(beta_i_c, na.rm = TRUE), .groups = "drop")

  chrom_depth_per_cell <- merge(
    chrom_depth_per_cell,
    beta_i_hat_per_chrom,
    by = c("chrom", "celltype"),
    all.x = TRUE
  )

  message("Step3: Estimating copy number factor")
  chrom_depth_per_cell$gamma_ci <-
    chrom_depth_per_cell$read_depth /
    (chrom_depth_per_cell$lib_size * chrom_depth_per_cell$beta_i_hat)

  message("Step4: Examining expected and observed depth")
  chrom_depth_per_cell$expected_depth <-
    chrom_depth_per_cell$lib_size * chrom_depth_per_cell$beta_i_hat
  chrom_depth_per_cell$standardized_residual <-
    (chrom_depth_per_cell$read_depth - chrom_depth_per_cell$expected_depth) /
    sqrt(chrom_depth_per_cell$expected_depth)

  chrom_depth_per_cell
}

#' Merge adjacent expected-depth bins with low cell counts
#'
#' Combines neighboring expected-depth bins until each merged bin contains more
#' than 200 observations, returning a mapping from original bins to merged bins.
#'
#' @param depth_bin_counts Data frame containing `depth_bin` and per-bin counts.
#'
#' @return A data frame mapping original expected-depth bins to merged bins.
merge_expected_depth_bins <- function(depth_bin_counts) {
  depth_bin_map <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(depth_bin_map) <- c("depth_bin", "depth_bin_new")
  n_rows <- nrow(depth_bin_counts)
  row_idx <- 1

  while (row_idx < n_rows) {
    if (depth_bin_counts[row_idx, 2] > 200) {
      depth_bin_map[nrow(depth_bin_map) + 1, ] <- c(depth_bin_counts[row_idx, ]$depth_bin[1], depth_bin_counts[row_idx, ]$depth_bin[1])
      row_idx <- row_idx + 1
    } else {
      for (next_idx in row_idx:nrow(depth_bin_counts)) {
        if (sum(depth_bin_counts[row_idx:next_idx, ]$n) > 200) {
          merged_rows <- as.data.frame(cbind(
            depth_bin_counts[row_idx:next_idx, ]$depth_bin,
            rep(depth_bin_counts[row_idx, ]$depth_bin[1], next_idx - row_idx + 1)
          ))
          colnames(merged_rows) <- c("depth_bin", "depth_bin_new")
          depth_bin_map <- rbind(depth_bin_map, merged_rows)
          row_idx <- next_idx + 1
          break
        }
        if (next_idx == nrow(depth_bin_counts)) {
          merged_rows <- as.data.frame(cbind(
            depth_bin_counts[row_idx:next_idx, ]$depth_bin,
            rep(depth_bin_counts[row_idx, ]$depth_bin[1], next_idx - row_idx + 1)
          ))
          colnames(merged_rows) <- c("depth_bin", "depth_bin_new")
          depth_bin_map <- rbind(depth_bin_map, merged_rows)
          row_idx <- next_idx + 1
          break
        }
      }
    }
  }

  depth_bin_map
}

#' Estimate a dispersion trend across mean-variance bins
#'
#' Fits a smooth trend to raw dispersion estimates and shrinks each bin-level
#' estimate toward the trend according to the number of observations in that bin.
#'
#' @param mean_variance_df Data frame containing per-bin means, variances, and
#'   observation counts.
#' @param prior_n Prior sample size controlling shrinkage toward the trend.
#'
#' @return The input table augmented with `phi_raw` and shrunken `phi`
#'   dispersion estimates.
estimate_dispersion_trend <- function(mean_variance_df, prior_n) {
  mean_variance_df$phi_raw <- (mean_variance_df$var - mean_variance_df$mean) / (mean_variance_df$mean^2)
  mean_variance_df$phi_raw <- 1 / mean_variance_df$phi_raw

  valid_rows <- is.finite(mean_variance_df$phi_raw) &
    mean_variance_df$phi_raw > 0 &
    is.finite(mean_variance_df$mean) &
    mean_variance_df$mean > 0
  if (sum(valid_rows) >= 5) {
    loess_fit <- stats::loess(
      log(phi_raw) ~ log(mean),
      data = mean_variance_df[valid_rows, ],
      span = 0.75,
      degree = 1
    )
    phi_trend <- exp(stats::predict(loess_fit, newdata = mean_variance_df))
  } else {
    phi_trend <- rep(stats::median(mean_variance_df$phi_raw[valid_rows], na.rm = TRUE), nrow(mean_variance_df))
  }

  shrinkage_weight <- mean_variance_df$n / (mean_variance_df$n + prior_n)
  mean_variance_df$phi <- exp((1 - shrinkage_weight) * log(phi_trend) + shrinkage_weight * log(mean_variance_df$phi_raw))

  mean_variance_df
}

#' Estimate mean-variance dispersion for chromosome-level CN testing
#'
#' Groups observations by expected depth, merges sparse depth bins, summarizes
#' their mean-variance relationship, and estimates a shrunken dispersion value
#' for each merged bin.
#'
#' @param chrom_depth_per_cell Output from [estimate_param()].
#' @param coverage_bin_size Relative width used to define expected-depth bins.
#' @param var_cut_off Reserved variance cutoff parameter retained for API
#'   compatibility.
#' @param prior_n Prior sample size used in shrinkage toward the dispersion
#'   trend.
#' @param min_bins Minimum number of expected-depth bins to construct.
#'
#' @return The input table augmented with merged depth-bin assignments and
#'   dispersion estimates.
estimate_variance_dispersion <- function(
  chrom_depth_per_cell,
  coverage_bin_size = 0.01,
  var_cut_off = 0.01,
  prior_n = 50,
  min_bins = 20
) {
  message("Step5: Grouping cells into bins by expected depth")

  max_expected_depth <- max(chrom_depth_per_cell$expected_depth)
  expected_depth_bin_width <- coverage_bin_size * max_expected_depth
  chrom_depth_per_cell$depth_bin <- round(chrom_depth_per_cell$expected_depth / expected_depth_bin_width, 0)

  num_bins <- length(unique(chrom_depth_per_cell$depth_bin))
  if (num_bins < min_bins) {
    expected_depth_bin_width <- max_expected_depth / min_bins
    chrom_depth_per_cell$depth_bin <- round(chrom_depth_per_cell$expected_depth / expected_depth_bin_width, 0)
  }

  depth_bin_counts <- chrom_depth_per_cell %>%
    dplyr::group_by(depth_bin) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(depth_bin)

  message("Step6: Merging neighbor bins if the number of cells in bin is < 200")
  depth_bin_map <- merge_expected_depth_bins(depth_bin_counts)

  chrom_depth_per_cell <- merge(chrom_depth_per_cell, depth_bin_map, by = "depth_bin")

  mean_variance_df <- chrom_depth_per_cell %>%
    dplyr::group_by(depth_bin_new) %>%
    dplyr::summarise(
      n = dplyr::n(),
      var = stats::var(read_depth, na.rm = TRUE),
      mean = mean(read_depth, na.rm = TRUE),
      .groups = "drop"
    )

  message("Step7: Estimating dispersion")
  mean_variance_df <- estimate_dispersion_trend(mean_variance_df, prior_n)

  chrom_depth_per_cell <- merge(chrom_depth_per_cell, mean_variance_df, by = "depth_bin_new", all.x = TRUE)

  chrom_depth_per_cell
}

#' Perform chromosome-level negative binomial CN testing
#'
#' Uses expected depth and dispersion estimates to compute per-chromosome
#' negative binomial p-values and copy-number calls for each cell.
#'
#' @param chrom_depth_per_cell Output from [estimate_variance_dispersion()].
#'
#' @return A data frame containing adjusted p-values and chromosome-level CN
#'   calls per cell.
cn_test_nb <- function(chrom_depth_per_cell) {
  message("Step8: Performing Negative Binomial test")

  chrom_depth_per_cell$p_value <- stats::pnbinom(
    q = chrom_depth_per_cell$read_depth,
    mu = chrom_depth_per_cell$lib_size * chrom_depth_per_cell$beta_i_hat,
    size = chrom_depth_per_cell$phi
  )

  chrom_depth_per_cell$p_value <- ifelse(
    chrom_depth_per_cell$read_depth > chrom_depth_per_cell$lib_size * chrom_depth_per_cell$beta_i_hat,
    1 - chrom_depth_per_cell$p_value,
    chrom_depth_per_cell$p_value
  )

  chrom_depth_per_cell$neg_log10_p <- -log10(chrom_depth_per_cell$p_value)
  chrom_depth_per_cell$p_value_adj <- stats::p.adjust(chrom_depth_per_cell$p_value, method = "BH")
  chrom_depth_per_cell$called_cna <- ifelse(chrom_depth_per_cell$p_value_adj < 0.01, "YES", "NO")
  chrom_depth_per_cell$cn_state_binom <- ifelse(
    chrom_depth_per_cell$p_value_adj < 0.01,
    chrom_depth_per_cell$gamma_ci * 2,
    NA_real_
  )

  chrom_depth_per_cell %>%
    dplyr::select(ID, chrom, p_value_adj, called_cna, cn_state_binom)
}

#' Run the full chromosome-level copy-number calling workflow
#'
#' Estimates expected depth parameters, dispersion, and chromosome-level
#' negative binomial CN calls from an input read-depth matrix.
#'
#' @param read_depth Data frame containing chromosome-level read depths and the
#'   required metadata columns.
#' @param positions Character vector of chromosome-level count columns.
#' @param min_lib_size Minimum library size required for a cell to be retained.
#' @param max_lib_size Maximum library size allowed for a cell to be retained.
#' @param coverage_bin_size Relative width used to define expected-depth bins.
#' @param var_cut_off Reserved variance cutoff parameter retained for API
#'   compatibility.
#'
#' @return A data frame of chromosome-level CN calls.
assign_cn_state.chrom <- function(
  read_depth,
  positions,
  min_lib_size = 5000,
  max_lib_size = 60000,
  coverage_bin_size = 0.01,
  var_cut_off = 0.01
) {
  chrom_depth_per_cell <- estimate_param(
    read_depth,
    positions,
    min_lib_size = min_lib_size,
    max_lib_size = max_lib_size
  )
  chrom_depth_per_cell <- estimate_variance_dispersion(
    chrom_depth_per_cell,
    coverage_bin_size = coverage_bin_size,
    var_cut_off = var_cut_off
  )
  cn_calls <- cn_test_nb(chrom_depth_per_cell)

  cn_calls
}
