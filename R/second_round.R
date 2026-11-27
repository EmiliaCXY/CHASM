# Second round: hybrid matrix-completion background on supported candidate regions, then strict re-segmentation.
#
# Why: robust PCA absorbs part of any CNA that many cells share into its low-rank "normal" background, so carriers lose
# signal as prevalence rises. The second round rebuilds the background over candidate regions (recurrent calls from a first
# round) from cells that are confidently CNA-free there, moves each cell's robust-PCA background to that level (keeping robust
# PCA's per-bin fit), and segments again with strict settings.
#
# Steps (chasm_second_round()):
#   1. round 1  : segment the robust-PCA residual at alpha_round1 (default 0.005) -> candidate regions and seed cells
#   2. strict   : segment the same residual at alpha_strict (default 0.001, V2 CN) -> confident calls
#   3. candidates (find_candidates()) and support (support_candidates(): binomial seed rule + optional overlap fix)
#   4. background (hybrid_background()): masked SVD on bins outside padded candidates -> cell factors; per candidate bin,
#      regression of the input on the factors over confident negatives; per cell and supported candidate, robust PCA's
#      background is shifted by (mean prediction - mean robust-PCA background) over the padded candidate
#   5. round 2  : re-segment P - L at alpha_strict (V2); a companion alpha_round1 call set updates the masking;
#      repeat 4-5 for `iterations` (candidates and support stay fixed)
# Final calls = round 2 at alpha_strict.

# bins table (variable, chrom, start, end, idx) in genomic order from padded bin names
.second_round_bins <- function(bin_names, chromosomes) {
  real <- bin_names[!grepl("padding", bin_names)]
  v <- sub("^read_depth_sqrt_centered\\.", "", real)
  p <- stringr::str_split_fixed(v, "-|_", 4)
  b <- data.frame(variable = v, chrom = p[, 1], start = as.numeric(p[, 2]), end = as.numeric(p[, 3]) + 1,
                  padded_name = real, stringsAsFactors = FALSE)
  b <- b[b$chrom %in% chromosomes, , drop = FALSE]
  b <- b[order(match(b$chrom, chromosomes), b$start), , drop = FALSE]
  b$idx <- seq_len(nrow(b)); rownames(b) <- NULL
  b
}

.nonzero_calls <- function(calls, cells) {
  x <- data.table::as.data.table(calls)[get("cn") != 2 & get("ID") %in% cells]
  unique(x[, c("ID", "variable", "cn"), with = FALSE])
}

#' Find candidate CNA regions from recurrent calls
#'
#' Candidate regions are runs of consecutive bins in which at least \code{min_cells} cells carry a call in the same
#' direction (sub-chromosomal calls only: chromosomes with more than one segment in that cell). A run is trimmed to bins
#' with at least \code{peak_frac} of the run's peak count and must keep at least \code{min_bins} bins.
#'
#' @param calls Data frame with \code{ID}, \code{variable}, \code{segment_id}, \code{cn} (e.g. \code{merge_calls()} output
#'   with \code{cn = cn_state_final}).
#' @param bins Data frame with \code{variable}, \code{chrom}, \code{start}, \code{end}, \code{idx} in genomic order.
#' @param min_cells,peak_frac,min_bins Candidate rules.
#'
#' @return A list of candidates, each a list with \code{chrom}, \code{start}, \code{end}, \code{dir} (-1 / 1), \code{bins},
#'   \code{peak}, \code{seeds} (cells called in that direction in the candidate) and \code{cn} (median CN of the seeds).
#' @export
find_candidates <- function(calls, bins, min_cells = 3L, peak_frac = 0.5, min_bins = 2L) {
  fin <- data.table::as.data.table(calls)
  if (!"chrom" %in% names(fin)) fin[, chrom := sub("_.*$", "", variable)]
  bins <- data.table::as.data.table(bins)
  f <- fin[cn != 2]
  multi <- fin[, list(nseg = data.table::uniqueN(segment_id)), by = c("ID", "chrom")][nseg > 1]
  f <- f[multi, on = c("ID", "chrom"), nomatch = 0][, dir := sign(cn - 2)]
  cnt <- f[, list(n = data.table::uniqueN(ID)), by = c("variable", "dir")]
  out <- list()
  for (dd in c(-1, 1)) {
    b <- data.table::copy(bins)[, n := 0L]; m <- cnt[dir == dd]; b[m, n := i.n, on = "variable"]
    b[, hit := n >= min_cells]
    b[, run := data.table::rleid(chrom, hit)]
    for (rr in unique(b[hit == TRUE, run])) {
      x <- b[run == rr]; pk <- max(x$n); x <- x[n >= peak_frac * pk]
      if (nrow(x) < min_bins) next
      cand_bins <- b[chrom == x$chrom[1] & idx >= min(x$idx) & idx <= max(x$idx), variable]
      seeds <- unique(f[dir == dd & variable %in% cand_bins, ID])
      out[[length(out) + 1]] <- list(chrom = x$chrom[1], start = min(x$start), end = max(x$end), dir = dd, bins = cand_bins, peak = pk,
                                     seeds = seeds, cn = round(stats::median(f[dir == dd & variable %in% cand_bins & ID %in% seeds, cn])))
    }
  }
  out
}

#' Decide which candidates get the second-round background
#'
#' Binomial seed rule: a candidate is supported if its seed count exceeds what one candidate could get if seeds were spread
#' at random at the median rate, \code{seeds >= qbinom(1 - fwer / #candidates, n_cells, median seeds / n_cells) + 1}.
#' Optional overlap fix: a supported candidate is dropped if it overlaps a supported candidate of the opposite direction
#' on the same chromosome and is at least twice as long (a broad artifact created by absorption of a strong focal event).
#'
#' @param cands Output of [find_candidates()].
#' @param n_cells Number of cells.
#' @param fwer Family-wise error rate over candidates.
#' @param drop_opposite_artifacts Apply the overlap fix.
#'
#' @return A data frame, one row per candidate: \code{seeds}, \code{threshold}, \code{supported_binomial}, \code{supported}.
#' @export
support_candidates <- function(cands, n_cells, fwer = 0.05, drop_opposite_artifacts = TRUE) {
  if (!length(cands)) return(data.frame(cand = integer(), seeds = integer(), threshold = numeric(), supported_binomial = logical(), supported = logical()))
  seed_n <- vapply(cands, function(cd) length(cd$seeds), integer(1))
  thr <- stats::qbinom(1 - fwer / length(cands), n_cells, stats::median(seed_n) / n_cells) + 1
  sup_b <- seed_n >= thr
  sup <- sup_b
  if (drop_opposite_artifacts) {
    for (i in which(sup_b)) for (j in which(sup_b)) {
      if (i != j && cands[[i]]$chrom == cands[[j]]$chrom && cands[[i]]$dir != cands[[j]]$dir &&
          length(intersect(cands[[i]]$bins, cands[[j]]$bins)) && length(cands[[i]]$bins) >= 2 * length(cands[[j]]$bins)) sup[i] <- FALSE
    }
  }
  data.frame(cand = seq_along(cands), chrom = vapply(cands, `[[`, "", "chrom"), start = vapply(cands, `[[`, 0, "start"),
             end = vapply(cands, `[[`, 0, "end"), dir = vapply(cands, `[[`, 0, "dir"), nbins = vapply(cands, function(cd) length(cd$bins), integer(1)),
             seeds = seed_n, threshold = thr, supported_binomial = sup_b, supported = sup, stringsAsFactors = FALSE)
}

# padded bins (regression) and window (negativity check) of each candidate
.candidate_geometry <- function(cands, bins, pad) {
  lapply(cands, function(cd_k) {
    cb <- bins[bins$chrom == cd_k$chrom, , drop = FALSE]; cin <- which(cb$variable %in% cd_k$bins); len <- length(cin)
    list(pad = cb$variable[max(1, min(cin) - pad):min(nrow(cb), max(cin) + pad)],
         win = cb$variable[max(1, min(cin) - max(len, pad)):min(nrow(cb), max(cin) + max(len, pad))])
  })
}

#' Hybrid second-round background
#'
#' @param P,L1 Cells x bins matrices (real bins, named by bin) of the input (sqrt-centred depth) and robust PCA's background.
#' @param geom Candidate geometry (padded bins and windows), from the internal helper.
#' @param supported Logical, one per candidate.
#' @param mask_calls Calls whose bins are masked (filled from \code{L1}) before the SVD (\code{ID}, \code{variable}).
#' @param neg_calls Calls defining non-negative cells: a cell with a call inside a candidate's window is not a negative there.
#' @param pos_calls Calls used as positives (with CN dummy terms) when \code{use_positives = TRUE}.
#' @param rank SVD rank (capped at half the number of bins outside the padded candidates).
#' @param use_positives Also fit positives (benchmarked configuration); default negatives only.
#' @param seed Seed set before the randomized SVD.
#'
#' @return A list: \code{L} (cells x bins background), \code{rank} used, \code{factor_bins}, \code{negatives} (per supported candidate).
#' @keywords internal
hybrid_background <- function(P, L1, geom, supported, mask_calls, neg_calls, pos_calls = NULL, rank = 10L, use_positives = FALSE, seed = 1234L) {
  if (!requireNamespace("rsvd", quietly = TRUE)) stop("hybrid_background: package 'rsvd' is required.", call. = FALSE)
  cells <- rownames(P)
  if (!any(supported)) return(list(L = L1, rank = NA_integer_, factor_bins = NA_integer_, negatives = data.frame()))
  RB <- unique(unlist(lapply(geom, `[[`, "pad"))); NC <- setdiff(colnames(P), RB)
  K <- min(as.integer(rank), length(NC) %/% 2L)
  if (K < 1L) {
    warning("hybrid_background: too few bins outside the padded candidates (", length(NC), "); second round skipped.", call. = FALSE)
    return(list(L = L1, rank = 0L, factor_bins = length(NC), negatives = data.frame()))
  }
  # masked SVD (hard-impute) on bins outside RB, each cell's own called bins masked and filled from L1
  W <- matrix(TRUE, nrow(P), length(NC), dimnames = list(cells, NC))
  o <- mask_calls[mask_calls$variable %in% NC & mask_calls$ID %in% cells, , drop = FALSE]
  W[cbind(match(o$ID, cells), match(o$variable, NC))] <- FALSE
  Z <- P[, NC, drop = FALSE]; Z[!W] <- L1[, NC, drop = FALSE][!W]
  set.seed(seed)
  for (it in 1:100) {
    sv <- rsvd::rsvd(Z, k = K); Lr <- sv$u %*% (sv$d * t(sv$v))
    ch <- sqrt(mean((Z[!W] - Lr[!W])^2)); Z[!W] <- Lr[!W]
    if (ch < 1e-4) break
  }
  U <- cbind(1, sv$u %*% diag(sv$d, nrow = length(sv$d)))
  # per candidate bin: regression on confident negatives (and optionally positives)
  near <- matrix(FALSE, nrow(P), length(RB), dimnames = list(cells, RB))
  for (g in geom) near[intersect(unique(neg_calls$ID[neg_calls$variable %in% g$win]), cells), g$pad] <- TRUE
  pos_cn <- matrix(2, nrow(P), length(RB), dimnames = list(cells, RB))   # strict calls: never negatives; positives only if use_positives
  if (!is.null(pos_calls)) {
    o <- pos_calls[pos_calls$variable %in% RB & pos_calls$ID %in% cells, , drop = FALSE]
    pos_cn[cbind(match(o$ID, cells), match(o$variable, RB))] <- o$cn
  }
  pred <- matrix(NA_real_, nrow(P), length(RB), dimnames = list(cells, RB))
  for (j in seq_along(RB)) {
    y <- P[, RB[j]]; pc <- pos_cn[, j]; isp <- pc != 2; isn <- !near[, j] & !isp
    lv <- sort(unique(pc[isp])); D <- if (use_positives && length(lv)) outer(pc, lv, "==") * 1 else NULL
    fitc <- if (use_positives) isp | isn else isn
    X <- cbind(U, D)
    b <- stats::lm.fit(X[fitc, , drop = FALSE], y[fitc])$coefficients; b[is.na(b)] <- 0
    bg <- drop(U %*% b[seq_len(ncol(U))]); bg <- bg + stats::median(y[isn] - bg[isn])
    pred[, j] <- bg
  }
  # hybrid: shift robust PCA's background per cell over each supported candidate (widest first; narrower ones overwrite)
  sh <- matrix(0, nrow(P), length(RB), dimnames = list(cells, RB))
  ks <- which(supported); ks <- ks[order(-vapply(geom[ks], function(g) length(g$pad), integer(1)))]
  for (kk in ks) {
    g <- geom[[kk]]
    sh[, g$pad] <- rowMeans(pred[, g$pad, drop = FALSE]) - rowMeans(L1[, g$pad, drop = FALSE])
  }
  L <- L1; L[, RB] <- L1[, RB] + sh
  negs <- data.frame(cand = which(supported), negatives = vapply(geom[supported], function(g) sum(!near[, g$pad[1]]), integer(1)))
  list(L = L, rank = K, factor_bins = length(NC), negatives = negs)
}

#' CHASM second round (hybrid background + strict segmentation)
#'
#' Runs round 1, the strict round, candidate finding and support, and the hybrid second round (see the description at the
#' top of \code{R/second_round.R}). Inputs are the objects of the standard pipeline.
#'
#' @param chrom_depth_per_cell Output of [normalize_depth()].
#' @param wavelet Output of [wavelet_transform()].
#' @param rpca Output of [robust_pca()] on \code{wavelet$mat.wavelet.transform}.
#' @param nb_calls Output of [assign_cn_state.chrom()] (chromosome-level NB track).
#' @param round1_calls Optional precomputed round-1 calls (\code{merge_calls()} output); computed when \code{NULL}.
#' @param alpha_round1,alpha_strict CBS alpha for round 1 / the strict round and round 2.
#' @param round1_backtransform Back-transform for round 1 (and the round-2 companion call set): "V2" (default, consistent with the
#'   strict round and round 2) or "V0" (as benchmarked). In a 12-setting comparison V2 here lowered detection of one-copy losses
#'   slightly (mean -0.015 vs the saved V0 calls) without fewer false calls.
#' @param negatives "any" (default, as benchmarked): confident negatives are cells with no round-1 or strict call in a candidate's
#'   window; "strict": no strict call only (mean -0.017 detection in the same comparison).
#' @param use_positives Also fit strict calls as positives with CN dummy terms (the benchmarked scripts did; no effect on detection in
#'   the comparison, 12 of 12 settings tied). Default \code{FALSE}.
#' @param drop_opposite_artifacts Overlap fix in [support_candidates()].
#' @param fwer,min_cells,pad,rank,iterations Tuning parameters (defaults as in the revision).
#' @param seed Base seed. Each segmentation and each masked SVD sets its own seed derived from it (round 1: seed + 1, strict: seed + 2,
#'   iteration i: SVD seed + 100 i, round-2 calls seed + 10 i + 1 / + 2), so results do not depend on the order of steps or on the session's
#'   random state. Robust PCA is computed outside this function; pass \code{robust_pca(..., seed = )} for a fully reproducible run.
#' @param chromosomes Chromosomes used.
#'
#' @return A list with \code{final} (round-2 calls at \code{alpha_strict}; \code{merge_calls()} format), \code{round1},
#'   \code{strict} (strict round, no second round), \code{candidates} ([support_candidates()] table), and \code{info}.
#' @export
chasm_second_round <- function(chrom_depth_per_cell, wavelet, rpca, nb_calls, round1_calls = NULL,
                               alpha_round1 = 0.005, alpha_strict = 0.001, round1_backtransform = c("V2", "V0"),
                               negatives = c("any", "strict"), use_positives = FALSE, drop_opposite_artifacts = TRUE,
                               fwer = 0.05, min_cells = 3L, pad = 2L, rank = 10L, iterations = 2L, seed = 1234L,
                               chromosomes = paste0("chr", c(1:22, "X", "Y"))) {
  round1_backtransform <- match.arg(round1_backtransform); negatives <- match.arg(negatives)
  f <- function(M) { M <- inv_wavelet_transform(M, wavelet$chrom.informed.wavelet); colnames(M) <- wavelet$bin.info.padded; M }
  Ppad <- f(wavelet$mat.wavelet.transform); L1pad <- f(rpca$Expected_Normal)
  id_of <- gsub("[#+]", "-", rownames(Ppad))          # call tables use this ID form; padded matrices keep the barcodes
  # each segmentation gets its own seed (CBS draws random numbers), so results do not depend on which calls ran before
  call_round <- function(Spad, Lpad, alpha, bt, k) {
    set.seed(seed + k)
    seg <- segment_residuals(Spad, alpha = alpha)
    cn <- assign_cn_state(chrom_depth_per_cell, Lpad, seg, backtransform = bt)
    out <- merge_calls(cn, nb_calls)
    out$chrom <- sub("_.*$", "", out$variable); out$cn <- out$cn_state_final
    out
  }
  bins <- .second_round_bins(colnames(Ppad), chromosomes)
  if (is.null(round1_calls)) round1_calls <- call_round(Ppad - L1pad, L1pad, alpha_round1, round1_backtransform, 1L)
  if (!"cn" %in% names(round1_calls)) round1_calls$cn <- round1_calls$cn_state_final
  if (!"chrom" %in% names(round1_calls)) round1_calls$chrom <- sub("_.*$", "", round1_calls$variable)
  strict <- call_round(Ppad - L1pad, L1pad, alpha_strict, "V2", 2L)

  rows <- which(id_of %in% unique(round1_calls$ID)); cells <- id_of[rows]          # cells kept by round 1 (e.g. NB library filter)
  col_of <- stats::setNames(match(bins$padded_name, colnames(Ppad)), bins$variable)
  P <- Ppad[rows, col_of, drop = FALSE]; L1 <- L1pad[rows, col_of, drop = FALSE]
  rownames(P) <- rownames(L1) <- cells; colnames(P) <- colnames(L1) <- bins$variable
  to_pad <- function(M, base) { out <- base; out[rows, col_of[colnames(M)]] <- M; out }

  cands <- find_candidates(round1_calls, bins, min_cells = min_cells)
  sup <- support_candidates(cands, length(cells), fwer = fwer, drop_opposite_artifacts = drop_opposite_artifacts)
  geom <- .candidate_geometry(cands, bins, pad)

  nz <- function(x) as.data.frame(.nonzero_calls(x, cells))
  conf <- nz(strict); anyc <- unique(rbind(nz(round1_calls), conf))
  final <- strict; info <- list()
  if (length(cands) && any(sup$supported)) {
    for (iter in seq_len(iterations)) {
      neg <- if (negatives == "strict") conf else anyc
      hb <- hybrid_background(P, L1, geom, sup$supported, mask_calls = anyc, neg_calls = neg, pos_calls = conf,
                              rank = rank, use_positives = use_positives, seed = seed + 100L * iter)
      S <- P - hb$L
      Spad <- to_pad(S, Ppad - L1pad); Lpad <- to_pad(hb$L, L1pad)
      c_lo <- call_round(Spad, Lpad, alpha_round1, round1_backtransform, 10L * iter + 1L)
      final <- call_round(Spad, Lpad, alpha_strict, "V2", 10L * iter + 2L)
      conf <- nz(final); anyc <- unique(rbind(nz(c_lo), conf))
      info[[iter]] <- list(iteration = iter, rank = hb$rank, factor_bins = hb$factor_bins, negatives = hb$negatives)
    }
  }
  drop_cols <- function(x) x[, setdiff(names(x), c("chrom", "cn")), drop = FALSE]
  list(final = drop_cols(final), round1 = drop_cols(round1_calls), strict = drop_cols(strict), candidates = sup, info = info)
}
