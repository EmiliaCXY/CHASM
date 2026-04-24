#' Merge wavelet-based and negative-binomial copy-number calls
#'
#' Joins segment-level CN calls from the wavelet pipeline with chromosome-level
#' negative-binomial calls and derives a final rounded copy-number state.
#'
#' @param cn_wavelet Data frame returned by [assign_cn_state()].
#' @param cn_nb Data frame returned by [assign_cn_state.chrom()].
#'
#' @return A merged data frame containing both call sets and a final
#'   `cn_state_final` column.
merge_calls <- function(cn_wavelet, cn_nb) {
  cn_nb$ID <- sub("\\+", "-", cn_nb$ID)
  cn_merged <- merge(
    cn_wavelet,
    cn_nb,
    by.x = c("ID", "chrom"),
    by.y = c("ID", "chrom"),
    suffixes = c("_wl", "_nb")
  )

  cn_merged$cn_state_final <- ifelse(
    cn_merged$called_cna == "YES" & !is.na(cn_merged$p_value_adj_wl),
    cn_merged$cn_state_binom,
    ifelse(
      cn_merged$called_cna == "NO" & !is.na(cn_merged$p_value_adj_wl) & cn_merged$p_value_adj_wl < 0.05,
      2,
      cn_merged$cn_state_adj_wl
    )
  )
  cn_merged$cn_state_final <- round(cn_merged$cn_state_final, 0)

  cn_merged
}
