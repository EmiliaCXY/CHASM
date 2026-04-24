estimate_param <- function(read_depth, positions) {
  message("Step1: Estimating cell specific factor: library size")

  read_depth <- as.data.frame(read_depth)
  if (!all(c("celltype", "ID") %in% colnames(read_depth))) {
    stop("estimate_param: read_depth must contain 'celltype' and 'ID' columns.")
  }
  if (!all(positions %in% colnames(read_depth))) {
    stop("estimate_param: positions not found in read_depth columns.")
  }

  read_depth$lib_size <- rowSums(read_depth[, positions, drop = FALSE], na.rm = TRUE)
  read_depth <- read_depth[read_depth$lib_size >= 5000, , drop = FALSE]
  read_depth <- read_depth[read_depth$lib_size <= 60000, , drop = FALSE]

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
  chrom_depth_per_cell$diff_normed2 <-
    (chrom_depth_per_cell$read_depth - chrom_depth_per_cell$expected_depth) /
    sqrt(chrom_depth_per_cell$expected_depth)

  chrom_depth_per_cell
}

estimate_variance_dispersion <- function(
  chrom_depth_per_cell,
  coverage_bin_size = 0.01,
  var_cut_off = 0.01,
  prior_n = 50,
  min_bins = 20
) {
  message("Step5: Grouping cells into bins by expected depth")

  max_expected_depth <- max(chrom_depth_per_cell$expected_depth)
  read_depth_bin_size <- coverage_bin_size * max_expected_depth
  chrom_depth_per_cell$depth_bin <- round(chrom_depth_per_cell$expected_depth / read_depth_bin_size, 0)

  num_bins <- length(unique(chrom_depth_per_cell$depth_bin))
  if (num_bins < min_bins) {
    read_depth_bin_size <- max_expected_depth / min_bins
    chrom_depth_per_cell$depth_bin <- round(chrom_depth_per_cell$expected_depth / read_depth_bin_size, 0)
  }

  sum_bin <- chrom_depth_per_cell %>%
    dplyr::group_by(depth_bin) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(depth_bin)

  message("Step6: Merging neighbor bins if the number of cells in bin is < 200")
  new_bin_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(new_bin_df) <- c("depth_bin", "depth_bin_new")
  num_row <- nrow(sum_bin)
  cur_row <- 1

  while (cur_row < num_row) {
    if (sum_bin[cur_row, 2] > 200) {
      new_bin_df[nrow(new_bin_df) + 1, ] <- c(sum_bin[cur_row, ]$depth_bin[1], sum_bin[cur_row, ]$depth_bin[1])
      cur_row <- cur_row + 1
    } else {
      for (next_row in cur_row:nrow(sum_bin)) {
        if (sum(sum_bin[cur_row:next_row, ]$n) > 200) {
          new_inserts <- as.data.frame(cbind(
            sum_bin[cur_row:next_row, ]$depth_bin,
            rep(sum_bin[cur_row, ]$depth_bin[1], next_row - cur_row + 1)
          ))
          colnames(new_inserts) <- c("depth_bin", "depth_bin_new")
          new_bin_df <- rbind(new_bin_df, new_inserts)
          cur_row <- next_row + 1
          break
        }
        if (next_row == nrow(sum_bin)) {
          new_inserts <- as.data.frame(cbind(
            sum_bin[cur_row:next_row, ]$depth_bin,
            rep(sum_bin[cur_row, ]$depth_bin[1], next_row - cur_row + 1)
          ))
          colnames(new_inserts) <- c("depth_bin", "depth_bin_new")
          new_bin_df <- rbind(new_bin_df, new_inserts)
          cur_row <- next_row + 1
          break
        }
      }
    }
  }

  chrom_depth_per_cell <- merge(chrom_depth_per_cell, new_bin_df, by = "depth_bin")

  mean_variance_df <- chrom_depth_per_cell %>%
    dplyr::group_by(depth_bin_new) %>%
    dplyr::summarise(
      n = dplyr::n(),
      var = stats::var(read_depth, na.rm = TRUE),
      mean = mean(read_depth, na.rm = TRUE),
      .groups = "drop"
    )

  message("Step7: Estimating dispersion")

  mean_variance_df$phi_raw <- (mean_variance_df$var - mean_variance_df$mean) / (mean_variance_df$mean^2)
  mean_variance_df$phi_raw <- 1 / mean_variance_df$phi_raw

  valid_idx <- is.finite(mean_variance_df$phi_raw) &
    mean_variance_df$phi_raw > 0 &
    is.finite(mean_variance_df$mean) &
    mean_variance_df$mean > 0
  if (sum(valid_idx) >= 5) {
    loess_fit <- stats::loess(
      log(phi_raw) ~ log(mean),
      data = mean_variance_df[valid_idx, ],
      span = 0.75,
      degree = 1
    )
    phi_trend <- exp(stats::predict(loess_fit, newdata = mean_variance_df))
  } else {
    phi_trend <- rep(stats::median(mean_variance_df$phi_raw[valid_idx], na.rm = TRUE), nrow(mean_variance_df))
  }

  w <- mean_variance_df$n / (mean_variance_df$n + prior_n)
  mean_variance_df$phi <- exp((1 - w) * log(phi_trend) + w * log(mean_variance_df$phi_raw))

  chrom_depth_per_cell <- merge(chrom_depth_per_cell, mean_variance_df, by = "depth_bin_new", all.x = TRUE)

  chrom_depth_per_cell
}

cn_test_nb <- function(chrom_depth_per_cell) {
  message("Step8: Performing Negative Binomial test")

  chrom_depth_per_cell$p_val <- stats::pnbinom(
    q = chrom_depth_per_cell$read_depth,
    mu = chrom_depth_per_cell$lib_size * chrom_depth_per_cell$beta_i_hat,
    size = chrom_depth_per_cell$phi
  )

  chrom_depth_per_cell$p_val <- ifelse(
    chrom_depth_per_cell$read_depth > chrom_depth_per_cell$lib_size * chrom_depth_per_cell$beta_i_hat,
    1 - chrom_depth_per_cell$p_val,
    chrom_depth_per_cell$p_val
  )

  chrom_depth_per_cell$neg_log10_p <- -log10(chrom_depth_per_cell$p_val)
  chrom_depth_per_cell$p_val_adjust <- stats::p.adjust(chrom_depth_per_cell$p_val, method = "BH")
  chrom_depth_per_cell$calledCNA <- ifelse(chrom_depth_per_cell$p_val_adjust < 0.05, "YES", "NO")
  chrom_depth_per_cell$CN_state.binom <- ifelse(
    chrom_depth_per_cell$p_val_adjust < 0.05,
    chrom_depth_per_cell$gamma_ci * 2,
    NA_real_
  )

  chrom_depth_per_cell %>%
    dplyr::select(ID, chrom, p_val_adjust, calledCNA, CN_state.binom)
}

assign_cn_state.chrom <- function(read_depth, positions) {
  chrom_depth_per_cell <- estimate_param(read_depth, positions)
  chrom_depth_per_cell <- estimate_variance_dispersion(chrom_depth_per_cell, coverage_bin_size = 0.01)
  cn_state_df <- cn_test_nb(chrom_depth_per_cell)

  cn_state_df
}
