#' Generate a pileup from BAM files
#'
#' @param Gene a character of gene name
#' @param regionsFile a region info file named SCISSOR_gaf.txt
#' @param BAMfiles a full path of BAM files from _Samtools_
#' @param caseIDs sample IDs
#' @param outputdir a directory to save outputs
#' @return a pileup matrix, regions, and ranges of genomic positions
#' @import SCISSOR
#' @export

gen_pileup <- function(Gene, regionsFile, BAMfiles, caseIDs, outputdir) {

  if (!Gene %in% regionsFile$gene_name) {
    stop(Gene, " is not in gene_name of SCISSOR_gaf.txt")
  }

  regions <- as.character(regionsFile[match(Gene, regionsFile$gene_name), c("regions")])
  Ranges = SCISSOR::get_Ranges(Gene=Gene, regions=regions, outputType="part_intron")
  pileup = SCISSOR::read_BAM(BAMfiles=BAMfiles, caseIDs=caseIDs, regions=regions)

  save(pileup, regions, Ranges, file=paste0(outputdir, Gene, "_pileup_part_intron.RData"))
}


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


#' Plot genome alignment profiles
#'
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param plot TRUE/FALSE turns on/off the genome alignment profiles plot. Default is TRUE.
#' @return a matrix and a plot, or a matrix for the percentages of sample properties where plot is TRUE or FALSE, respectively.
#' @import tidyverse dplyr ggplot2
#' @export

plot_GAP = function(sampleInfo, plot=TRUE) {

  # Calculate the percentages
  PCTmat0 <- sampleInfo %>%
    mutate(Unaligned=PF_BASES-PF_ALIGNED_BASES) %>%
    select(Unaligned, RIBOSOMAL_BASES, CODING_BASES, UTR_BASES, INTRONIC_BASES, INTERGENIC_BASES)

  PCTmat <- PCTmat0 %>%
    mutate(SUM=apply(PCTmat0, 1, sum),
           PCT_Unaligned=(Unaligned/SUM)*100,
           PCT_Intergenic=(INTERGENIC_BASES/SUM)*100,
           PCT_Intronic=(INTRONIC_BASES/SUM)*100,
           PCT_Coding.UTR=((CODING_BASES+UTR_BASES)/SUM)*100) %>%
    select(PCT_Unaligned, PCT_Intergenic, PCT_Intronic, PCT_Coding.UTR)

  PCTdf <- data.frame(grp=rep(c("Unaligned", "Intergenic", "Intronic", "Coding+UTR"), each=nrow(PCTmat)),
                      PCT=unlist(as.vector(PCTmat)))

  if (plot) {
    PCTdf$grp <- factor(PCTdf$grp, levels=c("Unaligned", "Intergenic", "Intronic", "Coding+UTR"), ordered=TRUE)

    p <- ggplot2::ggplot(PCTdf, aes(x=grp, y=PCT, color=grp)) +
      scale_colour_brewer(palette="Dark2") +
      geom_violin(position=position_dodge(width=0.9)) +
      geom_point(position=position_jitterdodge(seed=12345, dodge.width=0.9), alpha=0.4, shape=16, size=1) +
      stat_summary(fun=mean, geom="crossbar", width=0.5, position=position_dodge(width=0.9)) +
      labs(title="",x="", y="") +
      scale_x_discrete(limits=c("Unaligned", "Intergenic", "Intronic", "Coding+UTR")) +
      guides(colour=guide_legend(title="",nrow=1,byrow=TRUE)) +
      theme(legend.position="bottom",
            panel.background = element_rect(fill="gray97"))
    # print(p)
  }

  return(list(PCTmat=PCTmat, plot=p))
}


#' Filter low expressed genes
#'
#' @param genelist a vector of gene names
#' @param TPM a gene expression counts matrix transformed by TPM
#' @param thru threshold. Default is 5.
#' @param pct percent. Default is 50.
#' @return a vector of filtered gene names
#' @export

filter_lowExpGenes = function(genelist, TPM, thru=5, pct=50) {

  TPM.proteincoding <- na.omit(TPM[match(genelist, rownames(TPM)), ])
  rows_to_keep <- apply(TPM.proteincoding, 1, function(row) {mean(row<thru) < pct/100})
  genelist2 <- rownames(TPM.proteincoding[rows_to_keep, , drop=FALSE])

  return(genelist2)
}


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
#' @import parallel
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
#' @import
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
      pileupList[[g]] <- get_pileupExon(g, pileupPath)
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
  normlist = norm_pileup.list(pileupPath, geneNames=NULL, rnum=rnum, method=method)

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
  mar <- list(2, 3, 2:3)
  array_data <- simplify2array(scale.log.normlist)

  if (margin %in% 1:3) {
    # Use positive values only; if all values are 0 then all stats are 0
    pos_mask <- array_data>0
    var.sum <- apply(array_data, mar[[margin]], function(x) sum(x^2))
    pos_data <- replace(array_data, !pos_mask, NA)  # replace non-positive values with NA
    var.mean <- apply(pos_data, mar[[margin]], mean, na.rm=TRUE)
    var.sd <- apply(pos_data, mar[[margin]], sd, na.rm=TRUE)
    var.median <- apply(pos_data, mar[[margin]], median, na.rm=TRUE)
    var.mad <- apply(pos_data, mar[[margin]], mad, na.rm=TRUE)
    CV <- ifelse(var.mean<1e-10 | var.sd<1e-10, 0, var.sd/var.mean) # to adjust NaN, Inf
    robustCV <- ifelse(var.median<1e-10 | var.mad<1e-10, 0, var.mad/var.median)

  } else {
    stop(margin, " is not an option for margin.")
  }

  # Assemble the metrics
  if (margin %in% 1:2) {
    metrics <- data.frame(
      mean=var.mean,
      sd=var.sd,
      CV=CV,
      median=var.median,
      mad=var.mad,
      robustCV=robustCV
    )

  } else if (margin==3) {
    metrics <- data.frame(
      convert_pivot.longer(var.mean, c("sample", "gene", "mean")),
      sd=convert_pivot.longer(var.sd, c("sample", "gene", "sd"))[, "sd"],
      CV=convert_pivot.longer(CV, c("sample", "gene", "CV"))[, "CV"],
      median=convert_pivot.longer(var.median, c("sample", "gene", "median"))[, "median"],
      mad=convert_pivot.longer(var.mad, c("sample", "gene", "mad"))[, "mad"],
      robustCV=convert_pivot.longer(robustCV, c("sample", "gene", "robustCV"))[, "robustCV"]
    )
  }

  return(metrics)
}


#' Plot gene body coverage
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param scale TRUE/FALSE returns the scaled/unscaled normalized transcript coverage. Default is TRUE.
#' @param stat 1 and 2 return median and mean normalized coverage curves per sample, respectively. Default is 1.
#' @param plot TRUE/FALSE turns on/off the normalized transcript coverage plot. Default is TRUE.
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @return a matrix and a plot, or a matrix for the gene body coverage where plot is TRUE or FALSE, respectively.
#' @import tidyverse dplyr RColorBrewer ggplot2 matrixStats
#' @export

plot_GBC = function(pileupPath, geneNames, rnum=100, method=1, scale=TRUE, stat=2, plot=TRUE, sampleInfo) {

  scale.log.normlist = scale_pileup.list(pileupPath, geneNames, rnum=rnum, method=method, scale=TRUE)
  scale.arr <- simplify2array(scale.log.normlist)

  if (stat==1) {
    # Stat 1: Median curve per sample
    scale.geom <- matrix(matrixStats::rowMedians(matrix(scale.arr, nrow=prod(dim(scale.arr)[1:2]), ncol=dim(scale.arr)[3])), nrow=dim(scale.arr)[1])

  } else if (stat==2) {
    # Stat 2: Mean curve per sample
    scale.geom <- rowMeans(scale.arr, dims=2)

  } else {
    stop(stat," is not an option for stat.")
  }

  rownames(scale.geom) <- c(seq_len(rnum))
  colnames(scale.geom) <- colnames(scale.log.normlist[[1]])

  lgd <- sampleInfo %>%
    mutate(RatioIntron=INTRONIC_BASES/CODING_BASES) %>%
    select(SampleID, RatioIntron)

  GBP <- convert_Chr2Numcol(convert_pivot.longer(scale.geom, c("region", "sample", "scale.geom")) %>%
                              mutate(sample=gsub("\\.", "-", sample)) %>%
                              select(region, sample, scale.geom) %>%
                              inner_join(lgd, by=c("sample"="SampleID")), 1)

  if (plot) {
    myPalette <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))

    p <- ggplot2::ggplot(GBP, aes(x=region, colour=RatioIntron, group=sample)) +
      geom_line(aes(y=scale.geom), alpha=1, show.legend=TRUE) +
      labs(title="", x="Gene body percentile (5'\u21923')", y=paste0(c("Median", "Mean")[stat]," scaled normalized coverage")) +
      scale_colour_gradientn(colours=myPalette(100), limits=c(min(GBP$RatioIntron), max(GBP$RatioIntron)), name="Ratio intron") +
      theme(legend.position="bottom",
            panel.background=element_rect(fill="gray97"),
            strip.text.x=element_text(size=12, color="black"),
            strip.text.y=element_text(size=12, color="black"),
            strip.background=element_rect(color="NA", fill="white", linewidth=1, linetype="solid"),
            plot.title=element_text(hjust=0.5, face="bold"))
    # print(p)
  }

  return(list(GBP=GBP, plot=p))
}


#' Plot gene body coverage with good quality samples
#'
#' @param stat 1 and 2 return median and mean normalized coverage curves per sample, respectively. Default is 1.
#' @param plot TRUE/FALSE turns on/off the normalized transcript coverage plot. Default is TRUE.
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param GBCresult results of the gene body coverage with all samples
#' @param auc.vec a vector with SQI per sample
#' @return a matrix and a plot, or a matrix for the gene body coverage where plot is TRUE or FALSE, respectively.
#' @import tidyverse dplyr RColorBrewer ggplot2
#' @export

plot_GBCg = function(stat=2, plot=TRUE, sampleInfo, GBCresult, auc.vec) {

  # Update with good quality samples and PD
  GBP <- GBCresult$GBP %>%
    filter(sample %in% as.vector(auc.vec[auc.vec$SQI=="Good", c("Sample")])$Sample) %>%
    inner_join(auc.vec %>% select(Sample, PD), by=c("sample"="Sample"))

  if (plot) {
    myPalette <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))

    pRI <- ggplot2::ggplot(GBP, aes(x=region, colour=RatioIntron, group=sample)) +
      geom_line(aes(y=scale.geom), alpha=1, show.legend=TRUE) +
      labs(title="", x="Gene body percentile (5'\u21923')", y=paste0(c("Median", "Mean")[stat]," scaled normalized coverage")) +
      scale_colour_gradientn(colours=myPalette(100), limits=c(min(GBCresult$GBP$RatioIntron), max(GBCresult$GBP$RatioIntron)), name="Ratio intron") + # the original range of legend
      theme(legend.position="bottom",
            panel.background=element_rect(fill="gray97"),
            strip.text.x=element_text(size=12, color="black"),
            strip.text.y=element_text(size=12, color="black"),
            strip.background=element_rect(color="NA", fill="white", linewidth=1, linetype="solid"),
            plot.title=element_text(hjust=0.5, face="bold"))

    pPD <- ggplot2::ggplot(GBP, aes(x=region, colour=PD, group=sample)) +
      geom_line(aes(y=scale.geom), alpha=1, show.legend=TRUE) +
      labs(title="", x="Gene body percentile (5'\u21923')", y=paste0(c("Median", "Mean")[stat]," scaled normalized coverage")) +
      scale_colour_gradientn(colours=myPalette(100), limits=c(min(GBP$PD), max(GBP$PD)), name="PD") +
      theme(legend.position="bottom",
            panel.background=element_rect(fill="gray97"),
            strip.text.x=element_text(size=12, color="black"),
            strip.text.y=element_text(size=12, color="black"),
            strip.background=element_rect(color="NA", fill="white", linewidth=1, linetype="solid"),
            plot.title=element_text(hjust=0.5, face="bold"))

    # print(p)
  }

  return(list(GBP=GBP, plotRI=pRI, plotPD=pPD))
}


## -------------------------------------------------
## Sample Quality Index
## -------------------------------------------------

#' Get a mean coverage depth (MCD) for genes and samples
#'
#' @param genelist a vector of gene names
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @return MCD is a the number of genes x the number of samples matrix.
#' @import foreach doParallel SCISSOR parallel
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
#' @import foreach doParallel SCISSOR zoo parallel
#' @export

get_wCV <- function(genelist, pileupPath, sampleInfo, rnum=100, method=1, winSize=20, egPct=10) {

  if (!(2<=winSize && winSize<=rnum)) {
    stop("The window size ", winSize, " should be in [2, ", rnum, "].")
  }

  cl <- makeCluster(parallel::detectCores()-1)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add=TRUE)

  wCV <- foreach(g=1:length(pileupPath), .combine=rbind, .packages=c("SCISSOR", "zoo", "parallel"), .export=c("get_pileupExon", "norm_pileup.gene", "norm_pileup.spl")) %dopar%
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
#' @param rstPct restricted percent (one-side) to restrict genes by log transformed MC. Default is 20.
#' @param obsPct span includes the percent of observations in each local regression. Default is 50.
#' @return a vector with SQI per sample; a coordinate matrix of smoothed data; and a range of MCD.
#' @import ggplot2 DescTools SCISSOR tidyverse dplyr
#' @export

get_SQI = function(MCD, wCV, rstPct=20, obsPct=50) {

  auc.coord <- na.omit(data.frame(Gene=rep(rownames(MCD), ncol(MCD)),
                                  Sample=rep(colnames(MCD), each=nrow(MCD)),
                                  MCD=as.vector(MCD),
                                  wCV=as.vector(wCV))) %>%
    mutate(xMCD=log10(MCD+1)) %>%
    arrange(Sample, xMCD) # sort x-points for AUC

  # LOESS regression
  p <- ggplot2::ggplot(auc.coord, aes(x=xMCD, y=wCV)) +
    geom_smooth(data=auc.coord, aes(group=Sample), method="loess", span=obsPct/100, se=FALSE)
  smoothData <- ggplot_build(p)$data[[1]]

  # Map back to original group
  group_mapping <- auc.coord %>%
    distinct(Sample)
  group_mapping$group_id <- as.numeric(factor(group_mapping$Sample))
  smoothData <- smoothData %>%
    mutate(Sample=group_mapping$Sample[group])

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


#' Plot sample quality index (SQI) outputs
#'
#' @param SQIresult outputs from get_SQI function
#' @return figures for the distribution of SQI by PD; and the relation of wCV and MCD.
#' @references https://jtr13.github.io/cc21fall2/raincloud-plot-101-density-plot-or-boxplotwhy-not-do-both.html
#' @import ggplot2 tidyverse dplyr scales tibble ggpubr
#' @export

plot_SQI = function(SQIresult) {

  auc.vec2 <- result$auc.vec[, c("AUC", "PD")]
  auc.coord <- result$auc.coord
  rangeMCD <- result$rangeMCD
  df <- data.frame(var=rep(colnames(auc.vec2), each=nrow(auc.vec2)),
                   value=as.vector(as.matrix(auc.vec2)))
  df$var <- factor(df$var, levels=rev(colnames(auc.vec2)), ordered=TRUE)

  # Distribution of SQI by PD
  d <- df %>%
    ggplot2::ggplot(aes(x=var, y=value, fill=var)) +
    geom_flat_violin(position=position_nudge(x=0.2), alpha=0.4) +
    geom_point(aes(color=case_when(
      var=="AUC" ~ "grey",
      var=="PD" & value>3 ~ "red",
      var=="PD" & value<=3 ~ "darkgreen")),
      position=position_jitter(w=0.15, seed=12345), size=1.5, alpha=0.5, show.legend=F) +
    geom_boxplot(width=0.25, outlier.shape=NA, alpha=0.5) +
    scale_fill_manual(values=c(rep("grey", 2))) +
    scale_color_identity() +
    labs(x="", y="", fill="", title="") +
    guides(fill=guide_legend(nrow=1, byrow=TRUE))+
    theme(legend.position="none",
          legend.margin=margin(-5, 0, 0, 0),
          panel.background=element_rect(fill="gray97"),
          axis.text=element_text(size=25),
          axis.title=element_text(size=25)) +
    coord_flip()

  # Relation of wCV and MCD
  p <- ggplot2::ggplot(auc.coord, aes(x=x, y=y, group=Sample, color=SQI)) +
    geom_line(size=1) +
    scale_color_manual(values=c("Bad"=scales::alpha("red", 0.5), "Good"=scales::alpha("darkgreen", 0.5))) + # SQI per sample
    geom_rect(data=tibble::tibble(x1=rangeMCD[1], x2=rangeMCD[2], y1=-Inf, y2=+Inf),
              inherit.aes=FALSE,
              mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
              color="transparent",
              fill="blue",
              alpha=0.07) +
    xlab(expression(paste(log[10], "(MCD+1)"))) + ylab("wCV") +
    guides(size="none") +
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25),
          legend.title=element_text(size=18, face="bold"),
          legend.position=c(0.85, 0.85),
          legend.background=element_rect(colour=NA, fill=NA),
          legend.text=element_text(size=18),
          panel.background=element_rect(fill="gray97"))

  ggpubr::ggarrange(d, p, nrow=1, ncol=2, align="h")
}
