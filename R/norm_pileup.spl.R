#' Read depth normalization from pileup (for gene 1, sample 1)
#'
#' @param pileup a coverage pileup vector for one sample
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @return the normalized read depth is a vector with length=rnum.
#' @references https://github.com/hyochoi/SCISSOR
#' @import SCISSOR
#' @export

norm_pileup.spl = function(pileup, rnum=100, method=1) {

  row <- c(1:length(pileup))
  depthmat <- cbind(row, pileup)
  pos2 <- data.frame(round(seq(from=1, to=length(pileup), length.out=(2*rnum+1))))
  row_odd <- seq_len(nrow(pos2)) %% 2
  pos3 <- pos2[row_odd == 0, ] # even points
  pos4 <- pos2[row_odd == 1, ] # odd points
  region <- rep(1:rnum)

  if (method==1) {
    # Method 1: Raw value
    readdepth <- depthmat[pos3,2]

  } else if (method==2) {
    # Method 2: Interpolation
    readdepth <- 10^(rollmean(log10(depthmat[pos4,2]+1), 2))-1 # geometric mean

  } else {
    stop(method," is not an option for method.")
  }

  normmat <- as.matrix(cbind(region, pos3, readdepth))
  rownames(normmat) <- region
  colnames(normmat) <- c("region", "pos", "readdepth")

  return(normmat[,c("readdepth")])
}
