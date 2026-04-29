\name{plot_cn_heatmap}
\alias{plot_cn_heatmap}
\title{Plot a copy-number heatmap with chromosome gap separators}
\usage{
plot_cn_heatmap(
  cn_calls,
  cn_state_col = NULL,
  id_col = "ID",
  bin_col = "chrom_bin",
  chromosomes = c(paste0("chr", 1:22), "chrX", "chrY"),
  cluster_rows = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  max_cn = 4,
  min_cn = 0,
  color = grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
  file = NULL,
  ...
)
}
\arguments{
\item{cn_calls}{Data frame of per-bin copy-number calls, such as the output
of \code{\link[=merge_calls]{merge_calls()}}.}

\item{cn_state_col}{Name of the copy-number state column in \code{cn_calls}. If
\code{NULL}, the function uses \code{cn_state_final}, \code{cn_state_adj}, or
\code{CN_state.final} in that order.}

\item{id_col}{Name of the cell identifier column in \code{cn_calls}.}

\item{bin_col}{Name of the genomic bin column in \code{cn_calls}.}

\item{chromosomes}{Chromosome ordering to use when arranging bins.}

\item{cluster_rows}{Whether to cluster cells in the heatmap.}

\item{show_rownames}{Whether to display row names.}

\item{show_colnames}{Whether to display column names.}

\item{max_cn}{Maximum copy-number value to display in the heatmap.}

\item{min_cn}{Minimum copy-number value to display in the heatmap.}

\item{color}{Palette passed to \code{pheatmap::pheatmap()}.}

\item{file}{Optional file path passed to \code{pheatmap::pheatmap()}.}

\item{...}{Additional arguments passed to \code{pheatmap::pheatmap()}.}
}
\value{
A list containing the ordered matrix, chromosome gap positions, bin metadata,
and the \code{pheatmap} return object.
}
\description{
Reorders bins across chromosomes, reshapes merged copy-number calls into a
cell-by-bin matrix, and plots a heatmap with column gaps between
chromosomes.
}
