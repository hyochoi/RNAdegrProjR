#' Read depth normalization from pileupPath (for all genes, all samples)
#'
#' @param pileupPath file paths of coverage pileup .RData files that columns are samples for each file
#' @param geneNames file names of pileupData
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @return the normalized read depth is a rnum x the number of samples matrix at each gene list.
#' @references https://github.com/hyochoi/SCISSOR
#' @import SCISSOR miceadds
#' @export

norm_pileup.list = function(pileupPath, geneNames, rnum=100, method=1) {

  if(length(pileupPath)!=length(geneNames)) {
    stop("pileupPath must be the same length as geneNames")
  }

  if(length(pileupPath)==length(geneNames)) {
    pileupList <- list()
    for (g in 1:length(geneNames)){
      pileupList[[g]] <- list()
      PD <- miceadds::load.Rdata2(paste0(geneNames[g],".RData"), path=pileupPath[g])
      pileupList[[g]] <- PD
    }

    if (method==1 | method==2) {
      # Method 1: Raw value
      # Method 2: Interpolation
      normmat.list = lapply(pileupList, FUN=function(x) norm_pileup.gene(pileupData=x, rnum=rnum, method=method))

      # Add gene names in each list
      names(normmat.list) <- geneNames

    } else {
      stop(method," is not an option for method.")
    }
  }

  return(normmat.list)
}
