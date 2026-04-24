\name{assign_cn_state}
\alias{assign_cn_state}
\title{Assign copy-number states from denoised bin counts}
\description{Combines observed and expected counts with segmentation output to assign segment-level copy-number states.}
\usage{assign_cn_state(chrom_depth_per_cell, svd.normalized.read.depth, segment.output)}
\arguments{
  \item{chrom_depth_per_cell}{Long-format normalized bin data from \code{normalize_depth()}.}
  \item{svd.normalized.read.depth}{Denoised matrix in bin space.}
  \item{segment.output}{Segmentation output returned by \code{segment_residuals()}.}
}
\value{A data frame of segment-level copy-number state calls.}
