% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_normTC.R
\name{plot_normTC}
\alias{plot_normTC}
\title{Plot normalized transcript coverage using gene length normalization}
\usage{
plot_normTC(
  pileupPath,
  geneNames = NULL,
  rnum = 100,
  method = 1,
  scale = TRUE,
  stat = 1,
  plot = TRUE
)
}
\arguments{
\item{pileupPath}{file paths of coverage pileupData including .RData file names}

\item{geneNames}{gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.}

\item{rnum}{the number of regions for uniformly dividing the x-axis. Default is 100.}

\item{method}{1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.}

\item{scale}{TRUE/FALSE returns the scaled/unscaled normalized transcript coverage. Default is TRUE.}

\item{stat}{1 and 2 return median and mean normalized coverage curves per sample, respectively. Default is 1.}

\item{plot}{TRUE/FALSE turns on/off the normalized transcript coverage plot. Default is TRUE.}
}
\value{
a matrix and a plot, or a matrix for the normalized transcript coverage where plot is TRUE or FALSE, respectively.
}
\description{
Plot normalized transcript coverage using gene length normalization
}
\references{
https://github.com/hyochoi/SCISSOR
}
