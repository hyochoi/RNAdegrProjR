## -------------------------------------------------
##    code to review functions for gene length normalization
## -------------------------------------------------

## ---- renv

root <- "/Users/myeon/Documents/GitHub/RNAdegrProjR/"
devtools::load_all(root, quiet = TRUE)
devtools::load_all()

## ---- prepare the workspace

# imports
suppressPackageStartupMessages({
  library(dplyr) # For manipulating the data
  library(data.table) # For working with data.tables
  # library(biomaRt) # For getting gene information
  # library(org.Hs.eg.db) # For extracting entrezIDs
  library(SCISSOR) # For yaxis.hy
})


## ---- 1. pileData list for outputType="only_exon"
Allgenes <- list.files("data-raw/pileup_selected")
genes <- sub("_rs_pileup_part_intron.RData", "", Allgenes)
length(genes) # 101

num = c(1:length(genes))
path = sprintf("data-raw/pileupDataEx/PD%s.RData", num)


## ---- 2. Plots using functions (For KEAP1, all samples)
genes[84]
load(path[84])
pileupData <- PD84
AllPtLists <- colnames(pileupData)
length(AllPtLists) # 94


pdf(file="06_NormPileup_only_exon_output/norm_pileup_fnreview_KEAP1.pdf", width=12, height=9,
      bg="white", colormodel="cmyk", paper="letter")
  par(mfrow=c(4,2))

for (i in 1:length(AllPtLists)){
  pileup <- as.matrix(pileupData[,i])
  colnames(pileup) <- colnames(pileupData)[i]


  ## Make normmat
  row <- c(1:length(pileup))
  depthmat <- cbind(row, pileup)
  pos2 <- data.frame(round(seq(from=1, to=length(pileup), length.out=201)))
  row_odd <- seq_len(nrow(pos2)) %% 2
  pos3 <- pos2[row_odd==0, ] # even points
  pos4 <- pos2[row_odd==1, ] # odd points
  region <- rep(1:100)

  # Method 1: Raw value
  readdepth <- depthmat[pos3,2]

  normmat <- as.matrix(cbind(region, pos3, readdepth))
  rownames(normmat) <- region
  colnames(normmat) <- c("region", "pos", "readdepth")

  dim(normmat) # 100   3
  normmat1 <- normmat

  # Method 2: Interpolation
  readdepth <- 10^(rollmean(log10(depthmat[pos4,2]+1), 2))-1 # geometric mean

  normmat <- as.matrix(cbind(region, pos3, readdepth))
  rownames(normmat) <- region
  colnames(normmat) <- c("region", "pos", "readdepth")

  dim(normmat) # 100   3
  normmat2 <- normmat

  p1 <- substr(AllPtLists, start=1, stop=4)
  p2 <- substr(AllPtLists, start=6, stop=7)
  p3 <- substr(AllPtLists, start=9, stop=12)
  p4 <- substr(AllPtLists, start=14, stop=16)
  AllPtLists2 <- paste0(p1,".",p2,".",p3,".",p4) # length(AllPtLists2)=94


  ## Plot
  plot(pileupData[,AllPtLists[i]], type='l', lty=1, col="black", ylim=yaxis.hy(pileupData), main=sprintf('%s \nMethod 1: Raw value for KEAP1',AllPtLists[i]), ylab='', xlab='')
  abline(v=pos4,lty=3,col="red",lwd=1)
  points(x=normmat1[,2], y=normmat1[,3], col="green", pch=19, cex=0.5)

  plot(pileupData[,AllPtLists[i]], type='l', lty=1, col="black", ylim=yaxis.hy(pileupData), main=sprintf('%s \nMethod 2: Interpolation for KEAP1',AllPtLists[i]), ylab='', xlab='')
  abline(v=pos4,lty=3,col="red",lwd=1)
  points(x=pos4, y=depthmat[pos4,2], col="red", pch=19, cex=0.5)
  points(x=normmat2[,2], y=normmat2[,3], col="blue", pch=19, cex=0.5)


  mat1 = norm_pileup.spl(as.matrix(pileupData[,i]), 100, 1)
  plot(mat1, type='l', lty=1, col="green", ylim=yaxis.hy(pileupData), main=sprintf('%s \nnorm_pileup.spl for KEAP1',AllPtLists[i]), ylab='', xlab='')
  points(mat1, col="green", pch=19, cex=0.5)

  mat2 = norm_pileup.spl(as.matrix(pileupData[,i]), 100, 2)
  plot(mat2, type='l', lty=1, col="blue", ylim=yaxis.hy(pileupData), main=sprintf('%s \nnorm_pileup.spl for KEAP1',AllPtLists[i]), ylab='', xlab='')
  points(mat2, col="blue", pch=19, cex=0.5)


  normTC1 = plot_normTC(pileupPath=path, geneNames=genes, rnum=100, method=1, scale=TRUE, stat=1, plot=FALSE)
  normTC2 = normTC1[normTC1$sample==AllPtLists2[i], ]
  normTC3 = plot_normTC(pileupPath=path, geneNames=genes, rnum=100, method=2, scale=TRUE, stat=1, plot=FALSE)
  normTC4 = normTC3[normTC3$sample==AllPtLists2[i], ]
  normTC5 = plot_normTC(pileupPath=path, geneNames=genes, rnum=100, method=1, scale=TRUE, stat=2, plot=FALSE)
  normTC6 = normTC5[normTC5$sample==AllPtLists2[i], ]
  normTC7 = plot_normTC(pileupPath=path, geneNames=genes, rnum=100, method=2, scale=TRUE, stat=2, plot=FALSE)
  normTC8 = normTC7[normTC7$sample==AllPtLists2[i], ]

  plot(normTC2[,4], type='l', lty=1, col="green", ylim=yaxis.hy(normTC5[,4]), main=sprintf('%s \nscale.geomedian in plot_normTC for 101 selected genes',AllPtLists[i]), ylab='', xlab='')
  points(normTC2[,4], col="green", pch=19, cex=0.5)

  plot(normTC4[,4], type='l', lty=1, col="blue", ylim=yaxis.hy(normTC5[,4]), main=sprintf('%s \nscale.geomedian in plot_normTC for 101 selected genes',AllPtLists[i]), ylab='', xlab='')
  points(normTC4[,4], col="blue", pch=19, cex=0.5)

  plot(normTC6[,4], type='l', lty=1, col="green", ylim=yaxis.hy(normTC5[,4]), main=sprintf('%s \nscale.geomean in plot_normTC for 101 selected genes',AllPtLists[i]), ylab='', xlab='')
  points(normTC6[,4], col="green", pch=19, cex=0.5)

  plot(normTC8[,4], type='l', lty=1, col="blue", ylim=yaxis.hy(normTC5[,4]), main=sprintf('%s \nscale.geomean in plot_normTC for 101 selected genes',AllPtLists[i]), ylab='', xlab='')
  points(normTC8[,4], col="blue", pch=19, cex=0.5)


  print(i)
  }
  dev.off()
