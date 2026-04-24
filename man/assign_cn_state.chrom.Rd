\name{assign_cn_state.chrom}
\alias{assign_cn_state.chrom}
\title{Run the chromosome-level negative-binomial workflow}
\description{Convenience wrapper that estimates parameters, estimates dispersion, and returns chromosome-level CNA calls.}
\usage{assign_cn_state.chrom(read_depth, positions)}
\arguments{
  \item{read_depth}{Data frame with \code{ID}, \code{celltype}, and chromosome count columns.}
  \item{positions}{Character vector of chromosome columns to model.}
}
\value{A data frame of chromosome-level negative-binomial copy-number calls.}
