\name{segment_residuals}
\alias{segment_residuals}
\title{Segment residual signals with DNAcopy}
\description{Smooths per-bin residuals and segments them with \pkg{DNAcopy}.}
\usage{segment_residuals(residuals, alpha = 0.005)}
\arguments{
  \item{residuals}{Residual matrix with cells in rows and bins in columns.}
  \item{alpha}{Segmentation significance threshold passed to \code{DNAcopy::segment()}.}
}
\value{The \code{output} table from the segmented \code{DNAcopy} object.}
