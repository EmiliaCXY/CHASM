\name{construct_segment_to_bin_dictionary}
\alias{construct_segment_to_bin_dictionary}
\title{Map segmentation intervals back to bins}
\description{Uses genomic overlaps to connect segment intervals to the bins they span.}
\usage{construct_segment_to_bin_dictionary(output, bins)}
\arguments{
  \item{output}{Segmentation output with \code{chrom}, \code{loc.start}, \code{loc.end}, and \code{ID} columns.}
  \item{bins}{Character vector of genomic bin identifiers.}
}
\value{A data frame mapping segment identifiers to overlapping bins.}
