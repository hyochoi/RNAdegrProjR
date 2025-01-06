## -------------------------------------------------
## Comparison by RNA-seq Platforms
## -------------------------------------------------

#' Replace colnames of a matrix based on NewSampleId in Pilot2
#'
#' @param mat a matrix with entity_submitter_id columns
#' @param Pilot2 a reference matrix having NewSampleId
#' @export

rpl_NewSampleId = function(mat, Pilot2) {
  tmat <- t(mat)
  OldSampleId <- rownames(tmat)
  mat0 <- cbind(OldSampleId, tmat)

  # Keep the selected patients in Pilot data
  OldNewId <- Pilot2[, c("entity_submitter_id", "NewSampleId")]
  mat1 <- merge(OldNewId, mat0, by.x="entity_submitter_id", by.y="OldSampleId")

  # Change to a NewSampleId, sort by a NewSampleId
  mat2 <- mat1[order(mat1$NewSampleId), ]
  mat3 <- mat2[, 3:dim(mat2)[2]]
  rownames(mat3) <- c(mat2$NewSampleId)
  mat4 <- t(mat3)
  mat5 <- matrix(as.numeric(mat4), ncol=ncol(mat4))
  rownames(mat5) <- row.names(mat4)
  colnames(mat5) <- colnames(mat4)

  return(mat5)
}


#' Build coverage pileup based on annotations using NewSampleIds (modified build_pileup)
#'
#' @references https://github.com/hyochoi/SCISSOR
#' @import BiocManager Rsamtools SCISSOR
#' @export

build_pileupId = function(Pileup,caseIDs=NULL,regions,
                          inputType="part_intron",
                          outputType="part_intron") {
  # inputType  = "whole_intron", "part_intron", "only_exon"
  # output.ytpe = "whole_intron", "part_intron", "only_exon"

  if (missing(Pileup)) {
    stop("Pileup is missing")
  }
  if (missing(regions)) {
    stop("Regions should be specified")
  }

  strnd = strsplit(regions,":")[[1]][3]
  if (strnd=="-") {
    rawPileup = Pileup[rev(1:nrow(Pileup)),]
  } else {
    rawPileup = Pileup
  }
  if (is.null(caseIDs)) {
    #caseIDs = paste0("case-",1:ncol(rawPileup),sep="")
    caseIDs = colnames(Pileup)
  }
  if (inputType=="whole_intron") {

    if (outputType=="whole_intron") {
      intron.len = NULL
      ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
      covPileup = rawPileup;
    } else if (outputType=="part_intron") {
      intron.len = ceiling(len.intron.hy(regions)*0.5);
      ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
      covPileup = rawPileup[ep.new$coverage.col,] ;    #   Area to be included.
    } else if (outputType=="only_exon") {
      intron.len = 0;
      ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
      covPileup = rawPileup[ep.new$coverage.col,] ;    #   Area to be included.
    } else {
      stop(outputType," is not an option for outputType.")
    }

  } else if (inputType=="part_intron") {

    if (outputType=="whole_intron") {
      stop(outputType," is not an option when inputType=part_intron.")
    } else if (outputType=="part_intron") {
      intron.len = ceiling(len.intron.hy(regions)*0.5);
      ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
      covPileup = rawPileup ;    #   Area to be included.
      ##### If the dimension of the given pileup data is not equal to the expected one
      ##### from the intron.len calculated, should display an error. Fix this!
    } else if (outputType=="only_exon") {
      intron.len.temp = ceiling(len.intron.hy(regions)*0.5);
      ep.new.temp = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len.temp) ;

      exonic.region = c()
      for (i in 1:nrow(ep.new.temp$ep)) {
        exonic.region = c(exonic.region,(ep.new.temp$epl[i]:ep.new.temp$epr[i]));
      }
      covPileup = rawPileup[exonic.region,];
      rm(intron.len.temp,ep.new.temp,exonic.region);

      intron.len = 0;
      ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
    } else {
      stop(outputType," is not an option for outputType.")
    }

  } else if (inputType=="only_exon") {

    if (outputType=="whole_intron") {
      stop(outputType," is not an option when inputType==only_exon.")
    } else if (outputType=="part_intron") {
      stop(outputType," is not an option when inputType==only_exon.")
    } else if (outputType=="only_exon") {
      intron.len = 0;
      ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
      covPileup = rawPileup ;
    } else {
      stop(outputType," is not an option for outputType.")
    }

  } else {
    stop(inputType," is not an option for inputType.")
  }

  Ranges=get_Ranges(regions=regions,outputType=outputType)
  new.regions=Ranges$new.regions
  chr = strsplit(new.regions,":")[[1]][1]
  strtend = do.call(rbind,strsplit(strsplit(strsplit(new.regions,":")[[1]][2],",")[[1]],"-"))
  strnd = strsplit(new.regions,":")[[1]][3]
  strtend.num=matrix(as.numeric(strtend),ncol=2)
  allPos = unlist(sapply(1:nrow(strtend.num), function(x) strtend.num[x,1]:strtend.num[x,2]))

  rownames(covPileup) = allPos
  colnames(covPileup) = caseIDs
  if (strnd=="+") {
    return(covPileup)
  } else {
    return(covPileup[rev(1:nrow(covPileup)),])
  }
}


#' Plot pileup profiles with three pair lines (modified plot_pileup)
#'
#' @param AllPileup pileup for all samples to make the same y-axis with other samples
#' @param Pileup pileup for selected samples from AllPileup
#' @references https://github.com/hyochoi/SCISSOR
#' @import RColorBrewer SCISSOR
#' @export

plot_pileupPair = function(AllPileup,Pileup,Ranges,cases=NULL,logcount=NULL,
                           plot.meanpileup=TRUE,
                           col.pileup=NULL,col.meanpileup="grey",
                           main=NULL,cex.main=1.2,
                           print.ranges=TRUE,
                           xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,...) {

  ##  % needed variables
  if (missing(Pileup)) {
    stop("Pileup is missing")
  }
  if (missing(Ranges)) {
    stop("Genomic ranges should be needed (See get_Ranges)")
  }
  samplist = colnames(Pileup)
  n = ncol(Pileup)
  exons = matrix(Ranges$lRanges[,c(2,3)],ncol=2)
  if (is.null(cases)) cases = 1:n
  if (!is.numeric(cases)) {
    caseIDs=cases
    cases=which(samplist %in% caseIDs)
    if (length(cases)==0) {
      stop("No sample matches with the given case IDs.")
    }
  }

  candicol1 = c(brewer.pal(9,"Pastel1")[6], # candidate colors for exonic regions
                brewer.pal(8,"Pastel2")[6],
                brewer.pal(9,"YlOrBr")[1],
                brewer.pal(9,"YlOrRd")[1],
                brewer.pal(9,"YlOrRd")[2],
                brewer.pal(9,"Reds")[1],
                brewer.pal(9,"RdPu")[1],
                brewer.pal(9,"OrRd")[1],
                brewer.pal(9,"OrRd")[2],
                brewer.pal(9,"Oranges")[1],
                brewer.pal(9,"Oranges")[2],
                "aliceblue");
  candicol2 = brewer.pal(12,"Set3") # candidiate colors for regions with shape changes
  candicol3 = brewer.pal(8,"Pastel2") # candidiate colors for regions with shape changes
  candicol = c(candicol2,candicol3);
  exon.col = candicol1[9]

  # Set plot parameters
  if (is.null(xlim)) xlim = c(0,nrow(Pileup))
  #if (is.null(ylim)) ylim = yaxis.hy(Pileup)
  if (is.null(ylim)) ylim = yaxis.hy(AllPileup)
  if (is.null(xlab)) xlab = "Genomic positions"
  if (is.null(ylab)) ylab = "Read depth"
  if (is.null(main)) {
    if (length(cases)>1) {
      main = paste0(Ranges$Gene)
    } else if (length(cases)==1) {
      main = paste0(Ranges$Gene," | sample #",cases," (ID:",samplist[cases],")")
    }
  }
  # pileup colors
  # if (is.null(col.pileup)) {
  #   if (length(cases)>10) {
  #     x.mda = svd(Pileup) ;
  #     projmat = diag(x.mda$d)%*%t(x.mda$v) ;
  #     projmat[1,] = -projmat[1,] ;
  #     # colset = rainbow(n=length(projmat[1,]), start=0, end=0.756)[rank(-projmat[1,])]
  #     colset = colorRampPalette(brewer.pal(10, "Spectral"))(n)[rank(-projmat[1,])]
  #   } else {
  #     colset = rep("black",n)
  #   }
  # } else {
  #   if (length(col.pileup)==n) {
  #     colset=col.pileup
  #   } else {
  #     colset = rep("black",n)
  #     colset[cases] = rep(col.pileup,length=length(cases))
  #   }
  # }
  if (is.null(col.pileup)) {
    if (length(cases)==3) {
      colset = c("red","blue","green")
    }
  }

  # Start plotting
  # if (print.ranges) {
  #   par(mar=c(3.2,3.5,3,2))
  # } else {
  #   par(mar=c(3,3.5,3,2))
  # }
  par(mar=c(5,3.5,3,2))

  meanPileup = apply(Pileup, 1, median) ;
  plot(meanPileup, type='l', lty=2, lwd=0.5, ylim=c(min(0,ylim[1]),ylim[2]),
       xlim=xlim, axes=F, ylab=NA, xlab=NA, xaxs="i",yaxs="i", col="white") ;
  for (i in 1:nrow(exons)){
    #polygon(x=c(rep(exons[i,1],2),rep(exons[i,2],2)),y=c(-10000,(max(Pileup)+10000),(max(Pileup)+10000),-10000),col=exon.col,border=NA) ;
    polygon(x=c(rep(exons[i,1],2),rep(exons[i,2],2)),y=c(-10000,(max(AllPileup)+10000),(max(AllPileup)+10000),-10000),col=exon.col,border=NA) ;
  }
  abline(v=exons[,1],lty=1,col="lightyellow3",lwd=0.1) ;
  abline(v=exons[,2],lty=1,col="lightyellow3",lwd=0.1) ;
  title(main, cex.main=cex.main,font.main=1,line=0.5);

  # add legend
  legend("top", legend=c("FFM","FFT","PET"),
         col=c("red","blue","green"), lty=c(1, 1, 1), lwd=c(2,2,2), box.lty=0, cex=0.8,
         inset=c(0,1.2), xpd=TRUE, horiz=TRUE, x.intersp=0.5)

  if (print.ranges) {
    if (nrow(Ranges$lRanges)>1) {
      x.tick.at = c(1,Ranges$lRanges[2:nrow(Ranges$lRanges),2],max(Ranges$lRanges))
      x.labels.l = c(Ranges$Gene,Ranges$lRanges[2:nrow(Ranges$lRanges),2],max(Ranges$lRanges))
      x.labels.c = c(Ranges$Gene,Ranges$cRanges[2:nrow(Ranges$lRanges),1],max(Ranges$cRanges))
      x.labels.g = c(Ranges$chr,Ranges$gRanges[2:nrow(Ranges$gRanges),2],max(Ranges$gRanges))
      exon.tick.at = c(1,apply(Ranges$lRanges[,c(2,3)],1,mean))
    } else {
      x.tick.at = c(1,max(Ranges$lRanges))
      x.labels.l = c(Ranges$Gene,max(Ranges$lRanges))
      x.labels.c = c(Ranges$Gene,max(Ranges$cRanges))
      x.labels.g = c(Ranges$chr,max(Ranges$gRanges))
      exon.tick.at = c(1,(sum(Ranges$lRanges[,c(2,3)])*0.5))
    }
    exon.labels = c(Ranges$Gene,paste("E",1:dim(Ranges$lRanges)[1],sep=""))
    axis(side=1, tck=-0.01, at=x.tick.at, labels=NA, col.ticks="darkgrey") ;
    # axis(side=1, lwd=0, line=-1, cex.axis=0.8,col.axis="darkgrey",
    #      at=x.tick.at,labels=x.labels.c) ;
    axis(side=1, lwd=0, line=-1, cex.axis=0.8,col.axis="darkgrey",
         at=exon.tick.at,labels=exon.labels) ;
    axis(side=1, lwd=0, line=-0.1, cex.axis=0.8,col.axis="darkgrey",
         at=x.tick.at,labels=x.labels.g) ;
    mtext(side=1, xlab, line=2, cex=1) ;
  } else {
    mtext(side=1, xlab, line=1, cex=1)
  }
  if (!is.null(logcount)) {
    if (logcount==1) {
      labels = c(1,5,10,50,100,300,500,1000,2000,5000,10000,15000,20000,30000)
    } else {
      labels = c(5,10,50,100,300,500,1000,2000,5000,10000,15000,20000,30000)
    }
    tick.at = log10(labels+logcount)-log10(logcount);
  } else {
    tick.at = NULL;
    labels = TRUE;
  }
  axis(side=2, tck=-0.02, at=tick.at, col.ticks="darkgrey",las=1,
       labels=labels,lwd=0,line=-0.8,cex.axis=0.8,col.axis="darkgrey")
  mtext(side=2, ylab, line=2, cex=1) ;

  box(lwd=1.5)
  if (plot.meanpileup){
    points(meanPileup, type='l', lty=1, lwd=2, col=col.meanpileup) ;
  }
  for (case in cases) {
    points(Pileup[,case], type='l', lty=1, col=colset[case], ...) ;
  }
}


#' Plot median or mean pileup profiles within each group (modified plot_pileup)
#'
#' @param AllPileup pileup for all samples to make the same y-axis with other samples
#' @param Pileup pileup for selected samples from AllPileup
#' @param pair FFM or FFT or PET to set a different color line by pair
#' @references https://github.com/hyochoi/SCISSOR
#' @import RColorBrewer SCISSOR
#' @export

plot_pileupStat = function(AllPileup,Pileup,Ranges,cases=NULL,logcount=NULL,
                           plot.meanpileup=TRUE,
                           col.pileup=NULL,col.meanpileup="grey",
                           main=NULL,cex.main=1.2,
                           print.ranges=TRUE,
                           xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,pair,...) {

  ##  % needed variables
  if (missing(Pileup)) {
    stop("Pileup is missing")
  }
  if (missing(Ranges)) {
    stop("Genomic ranges should be needed (See get_Ranges)")
  }
  samplist = colnames(Pileup)
  n = ncol(Pileup)
  exons = matrix(Ranges$lRanges[,c(2,3)],ncol=2)
  # if (is.null(cases)) cases = 1:n
  # if (!is.numeric(cases)) {
  #   caseIDs=cases
  #   cases=which(samplist %in% caseIDs)
  #   if (length(cases)==0) {
  #     stop("No sample matches with the given case IDs.")
  #   }
  # }

  candicol1 = c(brewer.pal(9,"Pastel1")[6], # candidate colors for exonic regions
                brewer.pal(8,"Pastel2")[6],
                brewer.pal(9,"YlOrBr")[1],
                brewer.pal(9,"YlOrRd")[1],
                brewer.pal(9,"YlOrRd")[2],
                brewer.pal(9,"Reds")[1],
                brewer.pal(9,"RdPu")[1],
                brewer.pal(9,"OrRd")[1],
                brewer.pal(9,"OrRd")[2],
                brewer.pal(9,"Oranges")[1],
                brewer.pal(9,"Oranges")[2],
                "aliceblue");
  candicol2 = brewer.pal(12,"Set3") # candidiate colors for regions with shape changes
  candicol3 = brewer.pal(8,"Pastel2") # candidiate colors for regions with shape changes
  candicol = c(candicol2,candicol3);
  exon.col = candicol1[9]

  # Set plot parameters
  if (is.null(xlim)) xlim = c(0,nrow(Pileup))
  #if (is.null(ylim)) ylim = yaxis.hy(Pileup)
  if (is.null(ylim)) ylim = yaxis.hy(AllPileup)
  if (is.null(xlab)) xlab = "Genomic positions"
  if (is.null(ylab)) ylab = "Read depth"
  if (is.null(main)) {
    if (length(cases)>1) {
      main = paste0(Ranges$Gene)
    } else if (length(cases)==1) {
      main = paste0(Ranges$Gene," | sample #",cases," (ID:",samplist[cases],")")
    }
  }
  # pileup colors
  # if (is.null(col.pileup)) {
  #   if (length(cases)>10) {
  #     x.mda = svd(Pileup) ;
  #     projmat = diag(x.mda$d)%*%t(x.mda$v) ;
  #     projmat[1,] = -projmat[1,] ;
  #     # colset = rainbow(n=length(projmat[1,]), start=0, end=0.756)[rank(-projmat[1,])]
  #     colset = colorRampPalette(brewer.pal(10, "Spectral"))(n)[rank(-projmat[1,])]
  #   } else {
  #     colset = rep("black",n)
  #   }
  # } else {
  #   if (length(col.pileup)==n) {
  #     colset=col.pileup
  #   } else {
  #     colset = rep("black",n)
  #     colset[cases] = rep(col.pileup,length=length(cases))
  #   }
  # }

  if (pair=="All") {
    colset = c("red","blue","green")
  } else if (pair=="FFM") {
    colset = "red"
  } else if (pair=="FFT") {
    colset = "blue"
  } else if (pair=="PET") {
    colset = "green"
  } else {
    stop("The pair should be All or FFM or FFT or PET.")
  }

  # Start plotting
  # if (print.ranges) {
  #   par(mar=c(3.2,3.5,3,2))
  # } else {
  #   par(mar=c(3,3.5,3,2))
  # }
  par(mar=c(5,3.5,3,2))

  meanPileup = apply(Pileup, 1, median) ;
  plot(meanPileup, type='l', lty=2, lwd=0.5, ylim=c(min(0,ylim[1]),ylim[2]),
       xlim=xlim, axes=F, ylab=NA, xlab=NA, xaxs="i",yaxs="i", col="white") ;
  for (i in 1:nrow(exons)){
    #polygon(x=c(rep(exons[i,1],2),rep(exons[i,2],2)),y=c(-10000,(max(Pileup)+10000),(max(Pileup)+10000),-10000),col=exon.col,border=NA) ;
    polygon(x=c(rep(exons[i,1],2),rep(exons[i,2],2)),y=c(-10000,(max(AllPileup)+10000),(max(AllPileup)+10000),-10000),col=exon.col,border=NA) ;
  }
  abline(v=exons[,1],lty=1,col="lightyellow3",lwd=0.1) ;
  abline(v=exons[,2],lty=1,col="lightyellow3",lwd=0.1) ;
  title(main, cex.main=cex.main,font.main=1,line=0.5);

  # add legend
  legend("top", legend=c("FFM","FFT","PET"),
         col=c("red","blue","green"), lty=c(1, 1, 1), lwd=c(2,2,2), box.lty=0, cex=0.8,
         inset=c(0,1.2), xpd=TRUE, horiz=TRUE, x.intersp=0.5)

  if (print.ranges) {
    if (nrow(Ranges$lRanges)>1) {
      x.tick.at = c(1,Ranges$lRanges[2:nrow(Ranges$lRanges),2],max(Ranges$lRanges))
      x.labels.l = c(Ranges$Gene,Ranges$lRanges[2:nrow(Ranges$lRanges),2],max(Ranges$lRanges))
      x.labels.c = c(Ranges$Gene,Ranges$cRanges[2:nrow(Ranges$lRanges),1],max(Ranges$cRanges))
      x.labels.g = c(Ranges$chr,Ranges$gRanges[2:nrow(Ranges$gRanges),2],max(Ranges$gRanges))
      exon.tick.at = c(1,apply(Ranges$lRanges[,c(2,3)],1,mean))
    } else {
      x.tick.at = c(1,max(Ranges$lRanges))
      x.labels.l = c(Ranges$Gene,max(Ranges$lRanges))
      x.labels.c = c(Ranges$Gene,max(Ranges$cRanges))
      x.labels.g = c(Ranges$chr,max(Ranges$gRanges))
      exon.tick.at = c(1,(sum(Ranges$lRanges[,c(2,3)])*0.5))
    }
    exon.labels = c(Ranges$Gene,paste("E",1:dim(Ranges$lRanges)[1],sep=""))
    axis(side=1, tck=-0.01, at=x.tick.at, labels=NA, col.ticks="darkgrey") ;
    # axis(side=1, lwd=0, line=-1, cex.axis=0.8,col.axis="darkgrey",
    #      at=x.tick.at,labels=x.labels.c) ;
    axis(side=1, lwd=0, line=-1, cex.axis=0.8,col.axis="darkgrey",
         at=exon.tick.at,labels=exon.labels) ;
    axis(side=1, lwd=0, line=-0.1, cex.axis=0.8,col.axis="darkgrey",
         at=x.tick.at,labels=x.labels.g) ;
    mtext(side=1, xlab, line=2, cex=1) ;
  } else {
    mtext(side=1, xlab, line=1, cex=1)
  }
  if (!is.null(logcount)) {
    if (logcount==1) {
      labels = c(1,5,10,50,100,300,500,1000,2000,5000,10000,15000,20000,30000)
    } else {
      labels = c(5,10,50,100,300,500,1000,2000,5000,10000,15000,20000,30000)
    }
    tick.at = log10(labels+logcount)-log10(logcount);
  } else {
    tick.at = NULL;
    labels = TRUE;
  }
  axis(side=2, tck=-0.02, at=tick.at, col.ticks="darkgrey",las=1,
       labels=labels,lwd=0,line=-0.8,cex.axis=0.8,col.axis="darkgrey")
  mtext(side=2, ylab, line=2, cex=1) ;

  box(lwd=1.5)
  if (plot.meanpileup){
    points(meanPileup, type='l', lty=1, lwd=2, col=col.meanpileup) ;
  }

  for (case in cases) {
    #points(Pileup[,case], type='l', lty=1, col=colset[case], ...) ;
    points(Pileup[,case], type='l', lty=1, col=colset, ...) ;
  }

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
