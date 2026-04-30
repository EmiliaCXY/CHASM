\name{example_data_path}
\alias{example_data_path}
\title{Get the path to a bundled example data file}
\description{
Returns the path to one bundled example RDS file. During development, when
the package is sourced from the repository instead of installed, the function
falls back to the source-tree \code{inst/extdata} path.
}
\usage{
example_data_path(file_name)
}
\arguments{
  \item{file_name}{Name of the bundled example file. Use
  \code{\link[=list_example_data]{list_example_data()}} to see available options.}
}
\value{
A length-one character vector containing the path to the example RDS.
}
