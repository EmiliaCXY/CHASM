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

wavelet_transform <- function(chrom_depth_per_cell, bins, chromosomes) {
  if (!all(c("barcode", "bin", "read_depth_sqrt_centered") %in% colnames(chrom_depth_per_cell))) {
    stop("wavelet_transform: chrom_depth_per_cell must contain 'barcode', 'bin', 'read_depth_sqrt_centered'.")
  }

  message("Normalizing read depth and centering by expected bin factor...")

  read_depth_norm <- chrom_depth_per_cell[, c("barcode", "bin", "read_depth_sqrt_centered")]
  read_depth_norm_wide <- stats::reshape(
    read_depth_norm,
    idvar = "barcode",
    timevar = "bin",
    direction = "wide"
  )
  rownames(read_depth_norm_wide) <- read_depth_norm_wide$barcode
  read_depth_norm_wide <- subset(read_depth_norm_wide, select = -barcode)
  read_depth_norm_wide.mat <- as.matrix(read_depth_norm_wide)

  message("Reordering bins by chromosome and position...")

  bins_all <- colnames(read_depth_norm_wide)
  bins_df <- data.frame(
    bins = gsub("read_depth_sqrt_centered\\.", "", bins_all),
    stringsAsFactors = FALSE
  )
  bins_df[, c("chrom", "start", "end", "arm")] <- stringr::str_split_fixed(bins_df$bins, "-|_", 4)
  bins_df$chrom_arm <- paste0(bins_df$chrom, bins_df$arm)
  bins_df$start <- as.numeric(bins_df$start)

  message("Merging bin information with chromosome positions...")

  positions_df <- data.frame(chromosomes = chromosomes, order = seq_along(chromosomes))
  bins_df <- merge(positions_df, bins_df, by.x = "chromosomes", by.y = "chrom")
  bins_df <- bins_df %>% dplyr::arrange(order, chrom_arm, start)
  bins_df$bins <- paste0("read_depth_sqrt_centered.", bins_df$bins)

  read_depth_norm_wide.mat <- read_depth_norm_wide.mat[, bins_df$bins, drop = FALSE]

  message("Constructing chromosome-informed wavelet matrix...")

  bins_df$chrom_arm.adjust <- ifelse(
    bins_df$chromosomes %in% c("chr21", "chr22", "chrY"),
    bins_df$chromosomes,
    bins_df$chrom_arm
  )

  bins_num_per_chrom_arm <- bins_df %>%
    dplyr::group_by(chrom_arm.adjust, order) %>%
    dplyr::summarise(n.bins = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(order)

  bins_num_per_chrom_arm$n.bins.use <- 2 ^ ceiling(log2(bins_num_per_chrom_arm$n.bins))
  bins_num_per_chrom_arm <- bins_num_per_chrom_arm[bins_num_per_chrom_arm$n.bins.use > 2, ]

  total.bins <- sum(bins_num_per_chrom_arm$n.bins.use)
  chrom.informed.wavelet <- matrix(0, nrow = total.bins, ncol = total.bins)

  pre_row <- 1
  pre_col <- 1
  for (i in seq_len(nrow(bins_num_per_chrom_arm))) {
    chrom.bins.use <- bins_num_per_chrom_arm$n.bins.use[i]
    if (chrom.bins.use <= 2) {
      next
    }

    haarmat <- wavethresh::GenW(chrom.bins.use, filter.number = 1, family = "DaubExPhase")
    cur_row <- pre_row + chrom.bins.use - 1
    cur_col <- pre_col + chrom.bins.use - 1

    chrom.informed.wavelet[pre_row:cur_row, pre_col:cur_col] <- t(haarmat)

    pre_row <- cur_row + 1
    pre_col <- cur_col + 1
  }

  message("Padding read depth matrix to match wavelet dimensions...")

  read_depth_norm_wide.mat.new <- matrix(0, nrow = nrow(read_depth_norm_wide.mat), ncol = 0)
  rownames(read_depth_norm_wide.mat.new) <- rownames(read_depth_norm_wide.mat)

  for (i in seq_len(nrow(bins_num_per_chrom_arm))) {
    chrom_arm <- bins_num_per_chrom_arm$chrom_arm.adjust[i]
    chrom.bins.use <- bins_num_per_chrom_arm$n.bins.use[i]

    if (chrom.bins.use <= 2) {
      next
    }

    data.bins <- bins_df[bins_df$chrom_arm.adjust == chrom_arm, "bins"]
    data <- read_depth_norm_wide.mat[, data.bins, drop = FALSE]

    chrom.bins.padding <- chrom.bins.use - bins_num_per_chrom_arm$n.bins[i]
    if (chrom.bins.padding > 0) {
      padding <- matrix(0, nrow = nrow(read_depth_norm_wide.mat), ncol = chrom.bins.padding)
      colnames(padding) <- paste0(chrom_arm, ".padding.", seq_len(chrom.bins.padding))
      data.full <- cbind(data, padding)
    } else {
      data.full <- data
    }

    read_depth_norm_wide.mat.new <- cbind(read_depth_norm_wide.mat.new, data.full)
  }

  read_depth_norm_wide.mat.new.t <- t(read_depth_norm_wide.mat.new)

  message("Performing wavelet transform on read depth matrix...")

  mat.wavelet.transform <- chrom.informed.wavelet %*% read_depth_norm_wide.mat.new.t
  rownames(mat.wavelet.transform) <- rownames(read_depth_norm_wide.mat.new.t)

  list(
    mat.wavelet.transform = mat.wavelet.transform,
    chrom.informed.wavelet = chrom.informed.wavelet
  )
}

robust_pca <- function(mat.wavelet.transform) {
  if (!requireNamespace("rsvd", quietly = TRUE)) {
    stop("robust_pca: package 'rsvd' must be installed to run robust PCA.", call. = FALSE)
  }

  m <- mat.wavelet.transform
  lambda_val <- 1 / sqrt(max(nrow(m), ncol(m)))
  rpca_result <- rsvd::rrpca(m, lambda = lambda_val)

  l_hat <- rpca_result$L
  s_hat <- rpca_result$S

  colnames(l_hat) <- colnames(mat.wavelet.transform)
  colnames(s_hat) <- colnames(mat.wavelet.transform)

  list(
    Expected_Normal = l_hat,
    Sparse_Signal = s_hat
  )
}

inv_wavelet_transform <- function(mat.signal, mat.wavelet) {
  dat <- t(mat.wavelet) %*% mat.signal
  dat <- t(dat)

  rownames(dat) <- colnames(mat.signal)
  colnames(dat) <- rownames(mat.wavelet)

  dat
}
