#' Get metrics from scaled normalized transcript coverage
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param margin 1, 2, and 3 return metrics per sample, per gene, and across the genes per sample, respectively.
#' @return metrics including mean, sd, CV (sd/mean), median, mad, and robustCV (mad/median)
#' @references https://github.com/hyochoi/SCISSOR
#' @import SCISSOR tidyverse
#' @export

get_metrics = function(pileupPath, geneNames=NULL, rnum=100, method=1, margin) {

  scale.log.normarray = scale_pileup.list(path, genes, rnum=rnum, method=method)
  mar = list(2, 3, c(2,3))

  if (margin==1 | margin==2 | margin==3) {
    mean <- as.matrix(apply(scale.log.normarray, mar[[margin]], mean))
    sd <- as.matrix(apply(scale.log.normarray, mar[[margin]], sd))
    CV <- 0*mean
    CV[] <- ifelse((mean<1e-10 | sd<1e-10), 0, sd/mean) # to adjust NaN, Inf
    median <- as.matrix(apply(scale.log.normarray, mar[[margin]], median))
    mad <- as.matrix(apply(scale.log.normarray, mar[[margin]], mad))
    robustCV <- 0*median
    robustCV[] <- ifelse((median<1e-10 | mad<1e-10), 0, mad/median) # to adjust NaN, Inf

  } else {
    stop(margin," is not an option for margin.")
  }

  if (margin==1 | margin==2) {
    metrics <- data.frame(cbind(mean,sd,CV, median,mad,robustCV))
    colnames(metrics) <- c("mean","sd","CV", "median","mad","robustCV")

  } else if (margin==3) {
    vec.mean = convert_pivot.longer(mean, c("sample", "gene", "mean"))
    vec.sd = convert_pivot.longer(sd, c("sample", "gene", "sd"))
    vec.CV = convert_pivot.longer(CV, c("sample", "gene", "CV"))
    vec.median = convert_pivot.longer(median, c("sample", "gene", "median"))
    vec.mad = convert_pivot.longer(mad, c("sample", "gene", "mad"))
    vec.robustCV = convert_pivot.longer(robustCV, c("sample", "gene", "robustCV"))

    df <- data.frame(vec.mean, vec.sd[,c("sd")], vec.CV[,c("CV")], vec.median[,c("median")], vec.mad[,c("mad")], vec.robustCV[,c("robustCV")])
    metrics <- convert_Chr2Numcol(df, 3:8)
  }

  return(metrics)
}
