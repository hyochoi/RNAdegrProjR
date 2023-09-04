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
#' @export

get_metrics = function(pileupPath, geneNames=NULL, rnum=100, method=1, scale=TRUE, margin) {

  scale.log.normlist = scale_pileup.list(pileupPath, geneNames, rnum=rnum, method=method, scale=TRUE)
  mar = list(2, 3, 2:3)

  if (margin==1 | margin==2 | margin==3) {
    # Use positive values only; if all values are 0 then all stats are 0
    var.sum <- as.matrix(apply(simplify2array(scale.log.normlist), mar[[margin]], FUN=function(x) sum(x^2)))
    var.mean <- ifelse(var.sum==0, 0, as.matrix(apply(simplify2array(scale.log.normlist), mar[[margin]], FUN=function(x) mean(ifelse(x>0, x, NA), na.rm=TRUE))))
    var.sd <- ifelse(var.sum==0, 0, as.matrix(apply(simplify2array(scale.log.normlist), mar[[margin]], FUN=function(x) sd(ifelse(x>0, x, NA), na.rm=TRUE))))
    CV <- 0*var.mean
    CV[] <- ifelse((var.mean<1e-10 | var.sd<1e-10), 0, var.sd/var.mean) # to adjust NaN, Inf
    var.median <- ifelse(var.sum==0, 0, as.matrix(apply(simplify2array(scale.log.normlist), mar[[margin]], FUN=function(x) median(ifelse(x>0, x, NA), na.rm=TRUE))))
    var.mad <- ifelse(var.sum==0, 0, as.matrix(apply(simplify2array(scale.log.normlist), mar[[margin]], FUN=function(x) mad(ifelse(x>0, x, NA), na.rm=TRUE))))
    robustCV <- 0*var.median
    robustCV[] <- ifelse((var.median<1e-10 | var.mad<1e-10), 0, var.mad/var.median) # to adjust NaN, Inf

  } else {
    stop(margin," is not an option for margin.")
  }

  if (margin==1 | margin==2) {
    metrics <- data.frame(var.mean, var.sd, CV, var.median, var.mad, robustCV)
    colnames(metrics) <- c("mean", "sd", "CV", "median", "mad", "robustCV")

  } else if (margin==3) {
    vec.mean = convert_pivot.longer(var.mean, c("sample", "gene", "mean"))
    vec.sd = convert_pivot.longer(var.sd, c("sample", "gene", "sd"))
    vec.CV = convert_pivot.longer(CV, c("sample", "gene", "CV"))
    vec.median = convert_pivot.longer(var.median, c("sample", "gene", "median"))
    vec.mad = convert_pivot.longer(var.mad, c("sample", "gene", "mad"))
    vec.robustCV = convert_pivot.longer(robustCV, c("sample", "gene", "robustCV"))
    metrics <- data.frame(vec.mean, vec.sd[,c("sd")], vec.CV[,c("CV")], vec.median[,c("median")], vec.mad[,c("mad")], vec.robustCV[,c("robustCV")])
  }

  return(metrics)
}
