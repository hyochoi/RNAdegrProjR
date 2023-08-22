## -----------------------------------------
##    code for pileup plots for only exon
## -----------------------------------------

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
  library(SCISSOR) # For pileup and genomic ranges
})


## ---- 1. Pileup plots for each sample

load("data-raw/Pilot2.RData") # from DATASET.R

Allgenes <- list.files("data-raw/pileup_selected")
genes <- sub("_rs_pileup_part_intron.RData", "", Allgenes)
length(genes) # 101

for (g in 1:length(genes)){
  gthgene <- genes[g]
  load(file=sprintf("data-raw/pileup_selected/%s_rs_pileup_part_intron.RData", gthgene))
  pileup5 = rpl_NewSampleId(mat=pileup, Pilot2=Pilot2)

  ## Pileup figure 1 for all 101 genes (one pileup with three lines; only_exon; 34 samples)
  # Get genomic ranges for SCISSOR analysis
  # inputType  = "whole_intron", "part_intron", "only_exon"
  # outputType = "whole_intron", "part_intron", "only_exon"
  pileupData = build_pileupId(Pileup=pileup5, regions=regions, inputType="part_intron", outputType="only_exon")
  geneRanges = get_Ranges(Gene=Ranges$Gene, regions=regions, outputType="only_exon")

  # Plot for all patients
  AllPtLists <- substr(colnames(pileupData), start=1, stop=12)
  UnqPtLists <- unique(substr(colnames(pileupData), start=1, stop=12))

  pdf(file=sprintf("01_DataPreprocessing_output/gene_pileups_only_exon_34spls/%s.Fig1_plot_pileupPair.pdf", gthgene), width=12, height=9,
      bg="white", colormodel="cmyk", paper="letter")
  par(mfrow=c(3,3))
  for (i in 1:length(UnqPtLists)){
    plot_pileupPair(AllPileup=log10(pileupData+1), Pileup=log10(pileupData+1)[,(3*(i-1)+1):(3*(i-1)+3)], Ranges=geneRanges, main=AllPtLists[3*(i-1)+1], logcount=1)
  }
  dev.off()
}


## ---- 2. Median and mean pileup plots within each pair

for (g in 1:length(genes)){
  gthgene <- genes[g]
  load(file=sprintf("data-raw/pileup_selected/%s_rs_pileup_part_intron.RData", gthgene))
  pileup5 = rpl_NewSampleId(mat=pileup, Pilot2=Pilot2)

  # Remove 8 samples that had no expression in the pileup FFM (FF mRNA-seq)
  removespls <- c("TCGA-A6-2674-FFM","TCGA-A6-2674-FFT","TCGA-A6-2674-PET",
                  "TCGA-A6-2684-FFM","TCGA-A6-2684-FFT","TCGA-A6-2684-PET",
                  "TCGA-A6-3809-FFM","TCGA-A6-3809-FFT","TCGA-A6-3809-PET",
                  "TCGA-A6-3810-FFM","TCGA-A6-3810-FFT","TCGA-A6-3810-PET",
                  "TCGA-BK-A0CA-FFM","TCGA-BK-A0CA-FFT","TCGA-BK-A0CA-PET",
                  "TCGA-BK-A0CC-FFM","TCGA-BK-A0CC-FFT","TCGA-BK-A0CC-PET",
                  "TCGA-BK-A139-FFM","TCGA-BK-A139-FFT","TCGA-BK-A139-PET",
                  "TCGA-BK-A26L-FFM","TCGA-BK-A26L-FFT","TCGA-BK-A26L-PET")
  removecolnum = matrix(0, 1, 24)
  for (r in 1:24){
    removecolnum[,r] = which(grepl(pattern=removespls[r],x=colnames(pileup5)))
  }
  pileup6 <- pileup5[,-removecolnum]
  #dim(pileup5)[2] # 102
  #dim(pileup6)[2] # 78(102-8*3)

  # Get genomic ranges for SCISSOR analysis
  # inputType  = "whole_intron", "part_intron", "only_exon"
  # outputType = "whole_intron", "part_intron", "only_exon"
  pileupData = build_pileupId(Pileup=pileup6, regions=regions, inputType="part_intron", outputType="only_exon")
  geneRanges = get_Ranges(Gene=Ranges$Gene, regions=regions, outputType="only_exon")

  ## Pileup figure 2 for all 101 genes (median and mean pileup within each group; only_exon; 26 samples)
  FFMcol <- seq(1, length(AllPtLists), by=3)
  FFTcol <- FFMcol+1
  PETcol <- FFMcol+2

  pileupmFFM <- cbind(apply(pileupData[,FFMcol], 1, median), apply(pileupData[,FFMcol], 1, mean))
  pileupmFFT <- cbind(apply(pileupData[,FFTcol], 1, median), apply(pileupData[,FFTcol], 1, mean))
  pileupmPET <- cbind(apply(pileupData[,PETcol], 1, median), apply(pileupData[,PETcol], 1, mean))

  pdf(file=sprintf("01_DataPreprocessing_output/gene_pileups_only_exon_26spls/%s.Fig2_plot_pileupStat_only_exon_26spls.pdf", gthgene), width=12, height=9,
      bg="white", colormodel="cmyk", paper="letter")
  par(mfrow=c(3,3))
  plot_pileupStat(AllPileup=log10(pileupData+1), Pileup=log10(pileupmFFM+1), Ranges=geneRanges, case=1, main="FFM_median", pair="FFM", logcount=1)
  plot_pileupStat(AllPileup=log10(pileupData+1), Pileup=log10(pileupmFFT+1), Ranges=geneRanges, case=1, main="FFT_median", pair="FFT", logcount=1)
  plot_pileupStat(AllPileup=log10(pileupData+1), Pileup=log10(pileupmPET+1), Ranges=geneRanges, case=1, main="PET_median", pair="PET", logcount=1)

  plot_pileupStat(AllPileup=log10(pileupData+1), Pileup=log10(pileupmFFM+1), Ranges=geneRanges, case=2, main="FFM_mean", pair="FFM", logcount=1)
  plot_pileupStat(AllPileup=log10(pileupData+1), Pileup=log10(pileupmFFT+1), Ranges=geneRanges, case=2, main="FFT_mean", pair="FFT", logcount=1)
  plot_pileupStat(AllPileup=log10(pileupData+1), Pileup=log10(pileupmPET+1), Ranges=geneRanges, case=2, main="PET_mean", pair="PET", logcount=1)
  dev.off()
}
