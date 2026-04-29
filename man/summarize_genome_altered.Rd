\name{summarize_genome_altered}
\alias{summarize_genome_altered}
\title{Quantify the percentage of the genome with altered copy-number calls}
\usage{
summarize_genome_altered(
  cn_calls,
  bin_ids = NULL,
  id_col = "ID",
  segment_col = "segment_id",
  cn_state_col = NULL,
  diploid_state = 2
)
}
\arguments{
\item{cn_calls}{Data frame containing per-bin or repeated segment-level
copy-number assignments.}

\item{bin_ids}{Optional character vector of genomic bin identifiers. If
\code{NULL}, the function uses the \code{chrom_bin} column in \code{cn_calls}
when present, otherwise the \code{variable} column.}

\item{id_col}{Name of the cell identifier column in \code{cn_calls}.}

\item{segment_col}{Name of the segment identifier column in \code{cn_calls}.}

\item{cn_state_col}{Name of the copy-number state column in \code{cn_calls}. If
\code{NULL}, the function uses \code{cn_state_final}, \code{cn_state_adj}, or
\code{CN_state.final} in that order.}

\item{diploid_state}{Copy-number state treated as unaltered.}
}
\value{
A data frame with one row per cell and columns \code{ID},
\code{altered_length}, \code{genome_size}, and \code{pct_genome_altered}.
}
\description{
Summarizes per-cell copy-number output into the total altered genomic length
and percent of assayed bins assigned to a non-diploid copy-number state.
}
