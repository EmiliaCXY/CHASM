\name{recenter_svd_normalized_depth}
\alias{recenter_svd_normalized_depth}
\title{Recenter denoised bin values}
\description{Adds the expected bin-specific factor back to SVD- or RPCA-denoised bin values.}
\usage{recenter_svd_normalized_depth(svd.normalized.read.depth, chrom_depth_per_cell)}
\arguments{
  \item{svd.normalized.read.depth}{Denoised matrix in bin space.}
  \item{chrom_depth_per_cell}{Long-format output from \code{normalize_depth()}.}
}
\value{A recentered matrix aligned to the available bins.}
