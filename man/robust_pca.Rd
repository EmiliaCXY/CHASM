\name{robust_pca}
\alias{robust_pca}
\title{Run robust PCA on the wavelet matrix}
\description{Decomposes the wavelet-transformed matrix into low-rank expected signal and sparse residual signal.}
\usage{robust_pca(mat.wavelet.transform)}
\arguments{
  \item{mat.wavelet.transform}{Matrix returned by \code{wavelet_transform()}.}
}
\value{A list with \code{Expected_Normal} and \code{Sparse_Signal} matrices.}
\details{This function requires the optional \pkg{rsvd} package, which provides \code{rrpca()}.}
