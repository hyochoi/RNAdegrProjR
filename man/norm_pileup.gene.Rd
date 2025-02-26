% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm_pileup.gene.R
\name{norm_pileup.gene}
\alias{norm_pileup.gene}
\title{Gene length normalization for a pileup matrix (for gene 1, all samples)}
\usage{
norm_pileup.gene(pileupData, rnum = 100, method = 1)
}
\arguments{
\item{pileupData}{a coverage pileup matrix that columns are samples}

\item{rnum}{the number of regions for uniformly dividing the x-axis. Default is 100.}

\item{method}{1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.}
}
\value{
the normalized read depth is a rnum x the number of samples matrix.
}
\description{
Gene length normalization for a pileup matrix (for gene 1, all samples)
}
\references{
https://github.com/hyochoi/SCISSOR
}
