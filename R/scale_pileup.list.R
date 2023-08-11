#' Scale normalized transcript coverage using gene length normalization
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param scale TRUE/FALSE returns the scaled/unscaled normalized transcript coverage. Default is TRUE.
#' @return gene lists (rows are regions and columns are samples) for the normalized transcript coverage after gene length normalization
#' @references https://github.com/hyochoi/SCISSOR
#' @import SCISSOR
#' @export

scale_pileup.list = function(pileupPath, geneNames=NULL, rnum=100, method=1, scale=TRUE) {

  # Gene length normalization
  normlist = norm_pileup.list(path, genes, rnum=rnum, method=method)

  # Log-transformation
  log.normlist = lapply(normlist, FUN=function(x) log10(x+1))

  if (is.na(scale) | scale==TRUE) {
    scale.log.normlist = lapply(log.normlist, FUN=function(x) sweep(x, 2, apply(x,2,sum)+0.01, FUN="/"))
    return(scale.log.normlist)

  } else if (scale==FALSE) {
    return(log.normlist)

  } else {
    stop(scale," is not an option for scale.")
  }
}
