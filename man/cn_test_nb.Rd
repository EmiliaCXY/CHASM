\name{cn_test_nb}
\alias{cn_test_nb}
\title{Run chromosome-level negative-binomial testing}
\description{Computes one-sided negative-binomial tail probabilities, adjusts them for multiple testing, and derives chromosome-level CNA calls.}
\usage{cn_test_nb(chrom_depth_per_cell)}
\arguments{
  \item{chrom_depth_per_cell}{Output from \code{estimate_variance_dispersion()}.}
}
\value{A data frame of adjusted p-values and chromosome-level CNA calls.}
