#' Read depth normalization from pileupData (for gene 1, all samples)
#'
#' @param pileupData a coverage pileup matrix that columns are samples
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @return the normalized read depth is a rnum x the number of samples matrix.
#' @references https://github.com/hyochoi/SCISSOR
#' @import SCISSOR
#' @export

norm_pileup.gene = function(pileupData, rnum=100, method=1) {

  if (method==1 | method==2) {
    # Method 1: Raw value
    # Method 2: Interpolation
    normmat.gene = apply(pileupData, 2, FUN=function(x) norm_pileup.spl(pileup=x, rnum=rnum, method=method))

  } else {
    stop(method," is not an option for method.")
  }

  return(normmat.gene)
}
