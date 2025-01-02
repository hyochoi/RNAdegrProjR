## -------------------------------------------------
## Gene Length Normalization
## -------------------------------------------------

#' Gene length normalization for a single pileup (for gene 1, sample 1)
#'
#' @param pileup a coverage pileup vector for one sample
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @return the normalized read depth is a vector with length=rnum.
#' @references https://github.com/hyochoi/SCISSOR
#' @import zoo
#' @export

norm_pileup.spl = function(pileup, rnum=100, method=1) {

  row <- c(1:length(pileup))
  depthmat <- cbind(row, pileup)
  pos2 <- data.frame(round(seq(from=1, to=length(pileup), length.out=(2*rnum+1))))
  row_odd <- seq_len(nrow(pos2)) %% 2
  pos3 <- pos2[row_odd==0, ] # even points
  pos4 <- pos2[row_odd==1, ] # odd points
  region <- rep(1:rnum)

  if (method==1) {
    # Method 1: Raw value
    readdepth <- depthmat[pos3,2] # find read depth at even points (green points)

  } else if (method==2) {
    # Method 2: Interpolation
    readdepth <- 10^(round(zoo::rollmean(log10(depthmat[pos4,2]+1),2), 10))-1 # find geometric mean (blue points) using read depth at odd points (red points)

  } else {
    stop(method," is not an option for method.")
  }

  normmat <- as.matrix(cbind(region, pos3, readdepth))
  rownames(normmat) <- region
  colnames(normmat) <- c("region", "pos", "readdepth")

  return(normmat[,c("readdepth")])
}


#' Gene length normalization for a pileup matrix (for gene 1, all samples)
#'
#' @param pileupData a coverage pileup matrix that columns are samples
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @return the normalized read depth is a rnum x the number of samples matrix.
#' @references https://github.com/hyochoi/SCISSOR
#' @export

norm_pileup.gene = function(pileupData, rnum=100, method=1) {

  if (!(method %in% c(1, 2))) {
    stop(method, " is not an option for method.")
  }

  normmat.gene <- do.call(cbind, parallel::mclapply(seq_len(ncol(pileupData)), function(i) {
    norm_pileup.spl(pileup=pileupData[, i], rnum=rnum, method=method)
  }, mc.cores=parallel::detectCores()/2))
  colnames(normmat.gene) <- colnames(pileupData)

  return(normmat.gene)
}


#' Gene length normalization for pileup lists (for all genes, all samples)
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @return the normalized read depth is a rnum x the number of samples matrix at each gene list.
#' @references https://github.com/hyochoi/SCISSOR
#' @import miceadds
#' @export

norm_pileup.list = function(pileupPath, geneNames=NULL, rnum=100, method=1) {

  if (is.null(geneNames)) {
    geneNames = paste0("Gene_", c(1:length(pileupPath)))
  }

  if(length(pileupPath)!=length(geneNames)) {
    stop("pileupPath must be the same length as geneNames")

  } else if(length(pileupPath)==length(geneNames)) {
    pileupList <- list()
    for (g in 1:length(pileupPath)){
      pileupList[[g]] <- list()
      pileupList[[g]] <- miceadds::load.Rdata2(basename(pileupPath[g]), path=dirname(pileupPath[g]))
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


#' Scale normalized transcript coverage using gene length normalization
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param scale TRUE/FALSE returns the scaled/unscaled normalized transcript coverage. Default is TRUE.
#' @return gene lists (rows are regions and columns are samples) for the normalized transcript coverage after gene length normalization
#' @references https://github.com/hyochoi/SCISSOR
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
#' @import ggplot2
#' @export

plot_normTC = function(pileupPath, geneNames=NULL, rnum=100, method=1, scale=TRUE, stat=1, plot=TRUE) {

  scale.log.normlist = scale_pileup.list(pileupPath, geneNames, rnum=rnum, method=method, scale=TRUE)

  scale.geomedian <- as.matrix(apply(simplify2array(scale.log.normlist), 1:2, median)) # a rnum x the number of samples matrix
  scale.geomean <- as.matrix(apply(simplify2array(scale.log.normlist), 1:2, mean)) # a rnum x the number of samples matrix

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




## -------------------------------------------------
## Sample Quality Index
## -------------------------------------------------

#' Get a focused pileup of exon location for the g-th gene
#'
#' @param g the gene order in genelist
#' @param pileupPath file paths of coverage pileupData including .RData file names                                                     
#' @return a focused pileup is a the number of exon locations x the number of samples matrix for the g-th gene.
#' @import SCISSOR
#' @export

get_pileupExon = function(g, pileupPath) {

  load(file=pileupPath[g])

  # Keep exon location of union transcripts in pileup
  pileupData = SCISSOR::build_pileup(Pileup=pileup, regions=regions, inputType="part_intron", outputType="only_exon")
  colnames(pileupData) <- colnames(pileup) # to keep the original sample IDs

  return(pileupData)
}


#' Get a mean coverage depth (MCD) for genes and samples
#'
#' @param genelist a vector of gene names
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @return MCD is a the number of genes x the number of samples matrix.
#' @import SCISSOR
#' @export

get_MCD = function(genelist, pileupPath, sampleInfo) {
  cl <- makeCluster(parallel::detectCores()-1)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add=TRUE)

  MCD <- foreach(g=1:length(pileupPath), .combine=rbind, .packages=c("SCISSOR"), .export=c("get_pileupExon")) %dopar%
    {
      pileupData = get_pileupExon(g, pileupPath)
      if (nrow(pileupData) > 0) {
        # Use positive values only; if all values are 0 then all stats are 0
        sum <- colSums(pmax(pileupData, 0), na.rm=TRUE)
        n <- colSums(pileupData>0, na.rm=TRUE)
        ifelse(sum<1e-10 | n<1e-10, 0, sum/n) # to adjust NaN, Inf
      } else {
        rep(NA, nrow(sampleInfo))
      }
    }
  rownames(MCD) <- genelist

  return(MCD)
  stopCluster(cl)
}


#' Get a window coefficient of variation (wCV) for genes and samples
#'
#' @param genelist a vector of gene names
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param rnum the number of regions for uniformly dividing the x-axis for gene length normalization. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param winSize window size of the rolling window. Default is 20.
#' @param egPct edge percent (one-side) to calculate the trimmed mean. Default is 10.
#' @return wCV is a the number of genes x the number of samples matrix.
#' @import SCISSOR zoo
#' @export

get_wCV <- function(genelist, pileupPath, sampleInfo, rnum=100, method=1, winSize=20, egPct=10) {

  if (!(2<=winSize && winSize<=rnum)) {
    stop("The window size ", winSize, " should be in [2, ", rnum, "].")
  }

  cl <- makeCluster(parallel::detectCores()-1)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add=TRUE)

  wCV <- foreach(g=1:length(pileupPath), .combine=rbind, .packages = c("SCISSOR", "zoo", "parallel"), .export = c("get_pileupExon", "norm_pileup.gene", "norm_pileup.spl")) %dopar%
    {
      pileupData = get_pileupExon(g, pileupPath)
      if (nrow(pileupData) > 0) {
        # Normalization
        norm_pileup = norm_pileup.gene(pileupData, rnum=rnum, method=method)

        # Rolling CV in each sample
        rmean <- zoo::rollmean(norm_pileup, winSize, fill=NA, align="center", by.column=TRUE)
        rsd <- zoo::rollapply(norm_pileup, winSize, sd, fill=NA, align="center", by.column=TRUE)
        cv.mat <- ifelse(rmean<1e-10 | rsd<1e-10, 0, rsd/rmean) # to adjust NaN, Inf

        # 0-adjusted trimmed mean
        trmean_col <- function(column, trimFrac) {
          posVals <- column[column>0]
          ifelse(all(is.na(posVals)), NA, mean(posVals, na.rm=TRUE, trim=trimFrac))
        }

        wcv.vec <- unlist(parallel::mclapply(seq_len(ncol(cv.mat)), function(i) {
          trmean_col(cv.mat[, i], egPct/100)
        }, mc.cores=parallel::detectCores()/2))
        names(wcv.vec) <- colnames(cv.mat)
        wcv.vec

      } else {
        rep(NA, nrow(sampleInfo))
      }
    }
  rownames(wCV) <- genelist

  return(wCV)
  stopCluster(cl)
}

                                                     
#' Get a sample quality index (SQI) for samples
#'
#' @param MCD a mean coverage depth is a the number of genes x the number of samples matrix.
#' @param wCV a window coefficient of variation is a the number of genes x the number of samples matrix.
#' @param rstPct restricted percent (one-side) to restrict genes by log transformed MC. Default is 10.
#' @return a vector with SQI and a coordinate matrix by log transformed MCD, wCV, and smoothed wCV.
#' @import stats DescTools dplyr
#' @export

get_SQI = function(MCD, wCV, rstPct=20, obsPct=50) {

  auc.coord <- na.omit(data.frame(Gene=rep(rownames(MCD), ncol(MCD)),
                                  Sample=rep(colnames(MCD), each=nrow(MCD)),
                                  MCD=as.vector(MCD),
                                  wCV=as.vector(wCV))) %>%
    mutate(xMCD=log10(MCD+1)) %>%
    arrange(Sample, xMCD) # sort x-points for AUC

  # LOESS regression
  p <- ggplot(auc.coord, aes(x=xMCD, y=wCV)) +
    geom_smooth(data=auc.coord, aes(group=Sample), method="loess", span=obsPct/100, se=FALSE)
  smoothData <- ggplot_build(p)$data[[1]]

  # Map back to original group
  group_mapping <- auc.coord %>%
    distinct(Sample)
  group_mapping$group_id <- as.numeric(factor(group_mapping$Sample))
  smoothData <- smoothData %>%
    mutate(Sample = group_mapping$Sample[group])

  # Range of MCD
  posMCD <- MCD[MCD>0]
  rangeMin = log10(quantile(posMCD, probs=rstPct/100, na.rm=TRUE)+1)
  rangeMax = log10(quantile(posMCD, probs=1-rstPct/100, na.rm=TRUE)+1)

  auc.vec <- smoothData %>%
    filter(x>=rangeMin & x<rangeMax) %>% # restricted MCD
    group_by(Sample) %>%
    summarise(AUC=DescTools::AUC(x, y, method="spline")) %>% # calculate AUC
    mutate(PD=SCISSOR::pd.rate.hy(AUC, qrsc=TRUE), # projection depth
           SQI=ifelse(PD>3, "Bad", "Good")) # outlier detection

  auc.coord <- smoothData %>%
    select(x, y, Sample) %>%
    inner_join(auc.vec, by="Sample")

  return(list(auc.vec=auc.vec, auc.coord=auc.coord, rangeMCD=c(rangeMin, rangeMax)))
}
