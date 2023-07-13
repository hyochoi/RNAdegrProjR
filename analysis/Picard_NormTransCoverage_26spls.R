install.packages("ggstatsplot")
library(ggstatsplot)
install.packages("palmerpenguins")
library(palmerpenguins)
library(tidyverse)

install.packages("rstantools")
library(rstantools)
library(ggplot2)




###############################################################################
## Visualize the outputs from Picard tool CollectRnaSeqMetrics, HISTOGRAM
# Rows:
# Cols: normalized_position, NewSampleId
# Values: normalized_position, All_Reads.normalized_coverage
# Output: HSmat is a 101 x 103 matrix
###############################################################################
load("/Users/myeon/Documents/RNAdegr/output/01_DataPreprocessing_output/Pilot2.RData")

Allfiles <- list.files("/Users/myeon/Documents/RNAdegr/FFvsFFPE_ream_matrix")
files <- sub(".RNA_Metrics", "", Allfiles)

f=1
fthfile <- files[f]
all_content = readLines(sprintf("/Users/myeon/Documents/RNAdegr/FFvsFFPE_ream_matrix/%s.RNA_Metrics", fthfile))
keep_rows = all_content[c(12:112)]
histogram <- read.table(textConnection(keep_rows), sep="\t", header=F, na.strings="", stringsAsFactors=F)

Hmat = matrix(0, nrow(histogram),length(files)) # nrow(histogram)=101
for (f in 1:length(files)){
  fthfile <- files[f]
  all_content = readLines(sprintf("/Users/myeon/Documents/RNAdegr/FFvsFFPE_ream_matrix/%s.RNA_Metrics", fthfile))
  keep_rows = all_content[c(12:112)]
  histogram <- read.table(textConnection(keep_rows), sep="\t", header=F, na.strings="", stringsAsFactors=F)

  Hmat[,f] <- as.matrix(histogram[,2])
}
colnames(Hmat) <- c(files)

# Keep the selected patients in Pilot data
load_all()
Hmat1 = rpl_NewSampleId(mat=Hmat, Pilot2=Pilot2)

HSmat = cbind(histogram[,1], Hmat1)
colnames(HSmat) <- c("normalized_position", colnames(Hmat1))
dim(HSmat) # 101 103


##-----------------------------------------------------------------------
## Line plot with multiple lines
# Output: HSdf is a 303 x 28 data frame
#         Picard_NormTransCoverage_26spls.pdf
##-----------------------------------------------------------------------
# Remove 8 patients that had no expression in the pileup FFM (FF mRNA-seq)
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
  removecolnum[,r] = which(grepl(pattern=removespls[r],x=colnames(HSmat)))
}
HSmat0 <- HSmat[,-removecolnum]
dim(HSmat0) # 101  79[1+(102-8*3)]

FFMcol <- seq(2, dim(HSmat0)[2], by=3)
FFTcol <- FFMcol+1
PETcol <- FFMcol+2

pair2 <- c(rep("FFM",dim(HSmat0)[1]),rep("FFT",dim(HSmat0)[1]),rep("PET",dim(HSmat0)[1]))
FFMmat <- HSmat0[,c(1,FFMcol)]
FFTmat <- HSmat0[,c(1,FFTcol)]
PETmat <- HSmat0[,c(1,PETcol)]

pairmat <- rbind(FFMmat, FFTmat, PETmat)
HSdf <- data.frame(pair2, pairmat)
colnames(HSdf) <- c("pair2","normalized_position", substr(colnames(HSdf)[3:dim(HSdf)[2]], start=1, stop=12))
dim(HSdf) # 303(101*3)  28


# Line plot
ggplot(HSdf, aes(x=normalized_position, colour=pair2, group=pair2)) +
  geom_line(aes(y=HSdf[, 3 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 4 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 5 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 6 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 7 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 8 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 9 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 10 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 11 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 12 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 13 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 14 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 15 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 16 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 17 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 18 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 19 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 20 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 21 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 22 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 23 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 24 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 25 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 26 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 27 ]), alpha=0.4) +
  geom_line(aes(y=HSdf[, 28 ]), alpha=0.4) +
  geom_hline(yintercept=1, linetype=3) +
  scale_colour_brewer(palette="Set1") +
  theme(legend.position="bottom") +
  labs(title="",x="Normalized position", y="Mean normalized coverage") +
  guides(colour=guide_legend(title="")) +
  theme(panel.background = element_rect(fill="gray97"))

## end
