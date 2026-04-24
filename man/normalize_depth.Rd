\name{normalize_depth}
\alias{normalize_depth}
\title{Normalize per-bin read depth}
\description{Normalize per-bin read depth by library size and center it by an expected bin-specific factor.}
\usage{normalize_depth(df.depth, bins)}
\arguments{
  \item{df.depth}{Data frame with a \code{barcode} column and one column per genomic bin.}
  \item{bins}{Character vector of bin identifiers to normalize.}
}
\value{A long-format data frame with normalized and centered read-depth values per cell and bin.}
