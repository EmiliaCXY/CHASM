\name{inv_wavelet_transform}
\alias{inv_wavelet_transform}
\title{Invert the wavelet transformation}
\description{Projects a signal matrix from wavelet space back into bin space.}
\usage{inv_wavelet_transform(mat.signal, mat.wavelet)}
\arguments{
  \item{mat.signal}{Signal matrix in wavelet space.}
  \item{mat.wavelet}{Wavelet basis matrix returned by \code{wavelet_transform()}.}
}
\value{A matrix in bin space with cells in rows.}
