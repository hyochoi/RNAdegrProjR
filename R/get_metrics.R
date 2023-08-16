#' Get metrics from scaled normalized transcript coverage
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param scale TRUE/FALSE returns the scaled/unscaled normalized transcript coverage. Default is TRUE.
#' @param margin 1, 2, and 3 return metrics per sample, per gene, and across the genes per sample, respectively.
#' @return metrics including mean, sd, CV (sd/mean), median, mad, and robustCV (mad/median) per margin
#' @references https://github.com/hyochoi/SCISSOR
#' @import SCISSOR tidyverse
#' @export

get_metrics = function(pileupPath, geneNames=NULL, rnum=100, method=1, scale=TRUE, margin) {
  
  scale.log.normlist = scale_pileup.list(pileupPath, geneNames, rnum=rnum, method=method, scale=TRUE)
  mar = list(2, 3, 2:3)
  
  if (margin==1 | margin==2 | margin==3) {
    # Treat negative values as missing
    mean <- as.matrix(apply(simplify2array(scale.log.normlist), mar[[margin]], FUN=function(x) mean(ifelse((x>0 | x==0), x, NA), na.rm=TRUE)))
    sd <- as.matrix(apply(simplify2array(scale.log.normlist), mar[[margin]], FUN=function(x) sd(ifelse((x>0 | x==0), x, NA), na.rm=TRUE)))
    CV <- 0*mean
    CV[] <- ifelse((mean<1e-10 | sd<1e-10), 0, sd/mean) # to adjust NaN, Inf
    median <- as.matrix(apply(simplify2array(scale.log.normlist), mar[[margin]], FUN=function(x) median(ifelse((x>0 | x==0), x, NA), na.rm=TRUE)))
    mad <- as.matrix(apply(simplify2array(scale.log.normlist), mar[[margin]], FUN=function(x) mad(ifelse((x>0 | x==0), x, NA), na.rm=TRUE)))
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
    metrics <- data.frame(vec.mean, vec.sd[,c("sd")], vec.CV[,c("CV")], vec.median[,c("median")], vec.mad[,c("mad")], vec.robustCV[,c("robustCV")])
  }
  
  return(metrics)
}
