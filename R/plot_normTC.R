#' Plot normalized transcript coverage using gene length normalization
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param scale TRUE/FALSE returns the scaled/unscaled normalized transcript coverage. Default is TRUE.
#' @param stat 1 and 2 return median and mean normalized coverage curves per sample, respectively. Default is 1.
#' @param plot TRUE/FALSE turns on/off the normalized transcript coverage plot. Default is TRUE.
#' @return a matrix and a plot, or a matrix for the normalized transcript coverage where plot is TRUE or FALSE, respectively.
#' @references https://github.com/hyochoi/SCISSOR
#' @import SCISSOR tidyverse ggplot2
#' @export

plot_normTC = function(pileupPath, geneNames=NULL, rnum=100, method=1, scale=TRUE, stat=1, plot=TRUE) {
  
  scale.log.normlist = scale_pileup.list(pileupPath, geneNames, rnum=rnum, method=method, scale=TRUE)
  
  scale.geomedian <- as.matrix(apply(simplify2array(scale.log.normlist), 1:2, median))
  scale.geomean <- as.matrix(apply(simplify2array(scale.log.normlist), 1:2, mean))
  
  vec.scale.geomedian = convert_pivot.longer(scale.geomedian, c("region", "sample", "scale.geomedian"))
  vec.scale.geomean = convert_pivot.longer(scale.geomean, c("region", "sample", "scale.geomean"))
  
  pair <- substr(vec.scale.geomedian$sample, start=14, stop=16)
  df <- data.frame(pair, vec.scale.geomedian, vec.scale.geomean[,c("scale.geomean")])
  normTC <- convert_Chr2Numcol(df, 2)
  
  if (stat==1) {
    # Stat 1: Median curve per sample
    normTC1 <- normTC[,c("region", "sample", "pair", "scale.geomedian")]
    
  } else if (stat==2) {
    # Stat 2: Mean curve per sample
    normTC1 <- normTC[,c("region", "sample", "pair", "scale.geomean")]
    
  } else {
    stop(stat," is not an option for stat.")
  }
  
  if (plot) {
    mt <- c("Method 1: Raw value", "Method 2: Interpolation")
    st <- c("Median", "Mean")
    ntcplot <- ggplot2::ggplot(normTC1, aes(x=region, colour=pair, group=sample)) +
      geom_line(aes(y=normTC1[,4]), alpha=0.4) +
      scale_colour_brewer(palette="Set1") +
      theme(legend.position="bottom") +
      labs(title=mt[method],x="Regions", y=paste0(st[stat]," normalized coverage")) +
      guides(colour=guide_legend(title="")) +
      theme(panel.background = element_rect(fill="gray97"))
    print(ntcplot)
  }
  
  return(normTC1)
}
