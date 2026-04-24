\name{wavelet_transform}
\alias{wavelet_transform}
\title{Apply chromosome-informed wavelet transformation}
\description{Reorders bins by chromosome and arm, pads them to dyadic lengths, and applies a chromosome-informed Haar-like wavelet basis.}
\usage{wavelet_transform(chrom_depth_per_cell, bins, chromosomes)}
\arguments{
  \item{chrom_depth_per_cell}{Long-format output from \code{normalize_depth()}.}
  \item{bins}{Character vector of genomic bin identifiers.}
  \item{chromosomes}{Character vector specifying chromosome order.}
}
\value{A list with the transformed matrix and the wavelet basis matrix.}
