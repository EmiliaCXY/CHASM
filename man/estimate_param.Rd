\name{estimate_param}
\alias{estimate_param}
\title{Estimate chromosome-level negative-binomial parameters}
\description{Estimates library-size, expected chromosome fraction, copy-number ratio, and expected count terms for chromosome-level testing.}
\usage{estimate_param(read_depth, positions)}
\arguments{
  \item{read_depth}{Data frame with \code{ID}, \code{celltype}, and chromosome count columns.}
  \item{positions}{Character vector of chromosome columns to use.}
}
\value{A long-format data frame with estimated model parameters per cell and chromosome.}
