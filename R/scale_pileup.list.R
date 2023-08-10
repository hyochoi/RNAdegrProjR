#' Scale normalized transcript coverage using gene length normalization
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @return a scaled array (dimensions are regions, samples, genes) after gene length normalization
#' @references https://github.com/hyochoi/SCISSOR
#' @import SCISSOR abind
#' @export

scale_pileup.list = function(pileupPath, geneNames=NULL, rnum=100, method=1) {

  # Gene length normalization
  normlist = norm_pileup.list(path, genes, rnum=rnum, method=method)
  normarray <- abind::abind(normlist, along = 3)

  # Scale after log-transformation
  log.normarray = log10(normarray+1)
  colsum <- as.matrix(apply(log.normarray, c(2,3), sum))
  scale.log.normarray = 0*(log.normarray)
  for (i in 1:dim(normarray)[1]){
    for (j in 1:dim(normarray)[2]){
      for (k in 1:dim(normarray)[3]){
        scale.log.normarray[i,j,k] = log.normarray[i,j,k]/(colsum[j,k]+0.01)
      }
    }
  }

  return(scale.log.normarray)
}
