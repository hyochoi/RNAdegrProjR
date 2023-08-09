#' Plot normalized transcript coverage using gene length normalization
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param stat 1 and 2 return median and mean normalized coverage curves per sample, respectively. Default is 1.
#' @return a matrix and a plot for the normalized transcript coverage in list 1 and 2, respectively.
#' @references https://github.com/hyochoi/SCISSOR
#' @import SCISSOR abind ggplot2
#' @export

plot_normTC = function(pileupPath, geneNames=NULL, rnum=100, method=1, stat=1) {

  # Gene length normalization
  normlist = norm_pileup.list(path, genes, rnum=rnum, method=method)
  normarray <- abind::abind(normlist, along = 3)

  # Stat over gene list after log-transformation
  log.normarray = log10(normarray+1)
  geomedian <- as.matrix(apply(log.normarray, c(1,2), median))
  geomean <- as.matrix(apply(log.normarray, c(1,2), mean))
  vec.geomedian = as.vector(geomedian)
  vec.geomean = as.vector(geomean)

  # Scale
  colsum <- as.matrix(apply(log.normarray, c(2,3), sum))

  scale.log.normarray = 0*(log.normarray)
  for (i in 1:dim(normarray)[1]){
    for (j in 1:dim(normarray)[2]){
      for (k in 1:dim(normarray)[3]){
        scale.log.normarray[i,j,k] = log.normarray[i,j,k]/(colsum[j,k]+0.01)
      }
    }
  }
  scale.geomedian <- as.matrix(apply(scale.log.normarray, c(1,2), median))
  scale.geomean <- as.matrix(apply(scale.log.normarray, c(1,2), mean))
  vec.scale.geomedian = as.vector(scale.geomedian)
  vec.scale.geomean = as.vector(scale.geomean)

  region <- as.matrix(rep(c(1:rnum), dim(geomedian)[2]))
  sample <- repeach(as.matrix(colnames(geomedian)), rnum)
  pair <- substr(sample, start=14, stop=16)
  normTC <- data.frame(region, sample, pair, vec.geomedian, vec.geomean, vec.scale.geomedian, vec.scale.geomean)
  colnames(normTC) <- c("region", "sample", "pair", "geomedian", "geomean", "scale.geomedian", "scale.geomean")

  if (stat==1) {
    # Stat 1: Median curve per sample
    normTC1 <- normTC[,c("region", "sample", "pair", "scale.geomedian")]

  } else if (stat==2) {
    # Stat 2: Mean curve per sample
    normTC1 <- normTC[,c("region", "sample", "pair", "scale.geomean")]

  } else {
    stop(stat," is not an option for stat.")
  }

  mt <- c("Method 1: Raw value", "Method 2: Interpolation")
  st <- c("Median", "Mean")
  plot <- ggplot2::ggplot(normTC1, aes(x=region, colour=pair, group=sample)) +
    geom_line(aes(y=normTC1[,4]), alpha=0.4) +
    scale_colour_brewer(palette="Set1") +
    theme(legend.position="bottom") +
    labs(title=mt[method],x="Regions", y=paste0(st[stat]," normalized coverage")) +
    guides(colour=guide_legend(title="")) +
    theme(panel.background = element_rect(fill="gray97"))

  return(list(normTC1, plot))
}
