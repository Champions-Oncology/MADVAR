% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/madvar.R
\name{madvar}
\alias{madvar}
\title{Filter biological data by variance, or find a cutoff}
\usage{
madvar(data, mads = 2, must_genes = NULL, plot_density = FALSE, ...)
}
\arguments{
\item{data}{a numeric matrix or a numeric vector (with positive values). If matrix - Genes must be in rows}

\item{mads}{numeric (0-10]. Number of MADs to be added to the median for cutoff. Default=2 according to the outlier detection threshold.}

\item{must_genes}{a string vector with gene names that must be included in the output. Relevant only for matrix input.}

\item{plot_density}{logical. Whether to plot the density of data variance distribution with the calculated cutoffs.}

\item{...}{additional arguments to other functions}
}
\value{
If \code{plot_density} == FALSE: a filtered matrix (if \code{data} is a matrix) or a numeric cutoff (if \code{data} is a vector); \cr
If \code{plot_density} == TRUE: a \code{ggplot} object. The blue, gray and red vertical lines in the plot correspond to mode(Asselin estimate, bw=0.9), median, and final cutoff, respectively. \cr
}
\description{
This function calculates a cutoff based on data variance, aiming to filter out the "0 variance peak".
}
\examples{
\dontrun{
madvar(df, mads = 1, FALSE, must_genes = c("ERBB2", "BRCA1"), plot_density = TRUE)
}
}
