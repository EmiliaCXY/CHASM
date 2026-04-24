\name{estimate_variance_dispersion}
\alias{estimate_variance_dispersion}
\title{Estimate dispersion for chromosome-level counts}
\description{Bins chromosomes by expected count and estimates shrunken dispersion values for the negative-binomial model.}
\usage{
estimate_variance_dispersion(
  chrom_depth_per_cell,
  coverage_bin_size = 0.01,
  var_cut_off = 0.01,
  prior_n = 50,
  min_bins = 20
)
}
\arguments{
  \item{chrom_depth_per_cell}{Output from \code{estimate_param()}.}
  \item{coverage_bin_size}{Relative width of expected-depth bins.}
  \item{var_cut_off}{Retained for compatibility with the original script interface.}
  \item{prior_n}{Pseudo-count controlling shrinkage toward the mean-variance trend.}
  \item{min_bins}{Minimum number of expected-depth bins before neighbor merging.}
}
\value{The input data frame with variance and dispersion estimates appended.}
