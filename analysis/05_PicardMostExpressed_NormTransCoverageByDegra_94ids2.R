## ---------------------------------------------
##    code for normalized transcript coverage
## ---------------------------------------------

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
  library(RColorBrewer) # For continuous palettes
  library(ggplot2) # For ggplot
  library(gridExtra) # For multiple ggplots in pdf
})



#install.packages("ggstatsplot")
library(ggstatsplot)
#install.packages("palmerpenguins")
library(palmerpenguins)
library(tidyverse)

#install.packages("rstantools")
library(rstantools)


#install.packages("reshape2") # corr heatmap
library(reshape2)
library(devtools)



## ---- 1. Outputs from Picard tool CollectRnaSeqMetrics, HISTOGRAM

Allfiles <- list.files("data-raw/MostExpressed/1kb")
files <- sub("_1kb_5perc.RNA_Metrics", "", Allfiles)

f=1
fthfile <- files[f]
all_content = readLines(sprintf("data-raw/MostExpressed/1kb/%s_1kb_5perc.RNA_Metrics", fthfile))
keep_rows = all_content[c(12:112)]
histogram <- read.table(textConnection(keep_rows), sep="\t", header=F, na.strings="", stringsAsFactors=F)

# Make Hmat list from 4 length folders
AllLen <- c("1kb","1kb_5kb","5kb_10kb", "10kb")
L <- list()

for (l in 1:length(AllLen)){
  lthLen <- AllLen[l]
  L[[l]] <- list()

  Hmat = matrix(0, nrow(histogram),length(files)) # nrow(histogram)=101
  for (f in 1:length(files)){
    fthfile <- files[f]
    all_content = readLines(sprintf("data-raw/MostExpressed/%s/%s_%s_5perc.RNA_Metrics", lthLen,fthfile,lthLen))
    keep_rows = all_content[c(12:112)]
    histogram <- read.table(textConnection(keep_rows), sep="\t", header=F, na.strings="", stringsAsFactors=F)

    Hmat[,f] <- as.matrix(histogram[,2])
  }
  colnames(Hmat) <- c(files)
  L[[l]] <- Hmat
}
dim(L[[1]]) # 101  94
dim(L[[2]]) # 101  94
dim(L[[3]]) # 101  94
dim(L[[4]]) # 101  94

L1mat = cbind(histogram[,1], L[[1]])
L2mat = cbind(histogram[,1], L[[2]])
L3mat = cbind(histogram[,1], L[[3]])
L4mat = cbind(histogram[,1], L[[4]])
colnames(L1mat) <- c("normalized_position", colnames(L[[1]]))
colnames(L2mat) <- c("normalized_position", colnames(L[[2]]))
colnames(L3mat) <- c("normalized_position", colnames(L[[3]]))
colnames(L4mat) <- c("normalized_position", colnames(L[[4]]))
dim(L1mat) # 101  95
dim(L2mat) # 101  95
dim(L3mat) # 101  95
dim(L4mat) # 101  95

L1mat1 <- cbind(rep("1_0-1kb",nrow(L1mat)), L1mat)
L2mat1 <- cbind(rep("2_1-5kb",nrow(L2mat)), L2mat)
L3mat1 <- cbind(rep("3_5-10kb",nrow(L3mat)), L3mat)
L4mat1 <- cbind(rep("4_10+kb",nrow(L4mat)), L4mat)
HSmat <- rbind(L1mat1, L2mat1, L3mat1, L4mat1)
colnames(HSmat) <- c("GeneLen_kb", colnames(HSmat[,-1]))
dim(HSmat) # 404  96

# Separate by pair
pairspls <- c("FFM","FFT","PET")
paircol <- list()
pairmat <- list()
for (p in 1:length(pairspls)){
  paircol[[p]] <- list()
  pairmat[[p]] <- list()
  paircol[[p]] = grep(pattern=pairspls[p],x=colnames(HSmat))
  pairmat[[p]] <- HSmat[,c(1:2,paircol[[p]])]
}
dim(pairmat[[1]]) # 404  28
dim(pairmat[[2]]) # 404  36
dim(pairmat[[3]]) # 404  36


## ---- 2. Degradation variables from sample information table

load("data/SampleInfo.rda")

df <- list()
for (p in 1:length(pairspls)){
  df[[p]] <- list()

  pairmat1 <- pairmat[[p]]
  keepids = colnames(pairmat1[,3:ncol(pairmat1)])

  # Vectorize a matrix to a data frame
  len <- repmat(as.matrix(pairmat1[,1]), length(keepids), 1)
  normpos <- repmat(as.matrix(pairmat1[,2]), length(keepids), 1)
  ids <- repeach(as.matrix(keepids), nrow(pairmat1))
  pairdf1 <- cbind(len, normpos, ids, as.vector(pairmat1[,3:ncol(pairmat1)]))
  colnames(pairdf1) <- c("len", "normpos", "ids", "All_Reads.normalized_coverage")

  # Merge pairdf with degradation variables
  pairdf2 <- merge(pairdf1, SampleInfo[,c(1,37:39)], by.x="ids", by.y="NewSampleId", all.x=TRUE, sort=TRUE)
  pairdf3 = convert_Chr2Numcol(pairdf2,c(3:4))

  df[[p]] <- pairdf3
}
dim(df[[1]]) # 10504     7
dim(df[[2]]) # 13736     7
dim(df[[3]]) # 13736     7


## ---- 3. Normalized transcript coverage
# New facet label names
GeneLen_kb.labs <- c("0-1kb", "1-5kb", "5-10kb", "10+kb")
names(GeneLen_kb.labs) <- c("1_0-1kb", "2_1-5kb", "3_5-10kb", "4_10+kb")

# Continuous legends
myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))


## FFM where p=1
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(df[[1]]$ratio_28s_18s),max(df[[1]]$ratio_28s_18s)), name="Ratio 28s/18s")
Ratio28s18s1 <- ggplot(df[[1]], aes(x=normpos, colour=ratio_28s_18s, group=ids)) +
  geom_line(aes(y=df[[1]][, 4]), alpha=1, show.legend = TRUE) +
  geom_hline(yintercept=1, linetype=3) +
  labs(title=sprintf("%s",pairspls[1]),x="Normalized position", y="Mean normalized coverage") +
  theme(legend.position="bottom") +
  sc +
  theme(panel.background = element_rect(fill="gray97")) +
  coord_cartesian(ylim = c(0, 2)) +
  facet_wrap(~ len, ncol=2, labeller = labeller(len = GeneLen_kb.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid"))

sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(df[[1]]$rinvalue),max(df[[1]]$rinvalue)), name="RIN")
RIN1 <- ggplot(df[[1]], aes(x=normpos, colour=rinvalue, group=ids)) +
  geom_line(aes(y=df[[1]][, 4]), alpha=1, show.legend = TRUE) +
  geom_hline(yintercept=1, linetype=3) +
  labs(title=sprintf("%s",pairspls[1]),x="Normalized position", y="Mean normalized coverage") +
  theme(legend.position="bottom") +
  sc +
  theme(panel.background = element_rect(fill="gray97")) +
  coord_cartesian(ylim = c(0, 2)) +
  facet_wrap(~ len, ncol=2, labeller = labeller(len = GeneLen_kb.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid"))

sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(df[[1]]$RatioIntron),max(df[[1]]$RatioIntron)), name="Ratio INTRONIC BASES/CODING BASES")
RatioIntron1 <- ggplot(df[[1]], aes(x=normpos, colour=RatioIntron, group=ids)) +
  geom_line(aes(y=df[[1]][, 4]), alpha=1, show.legend = TRUE) +
  geom_hline(yintercept=1, linetype=3) +
  labs(title=sprintf("%s",pairspls[1]),x="Normalized position", y="Mean normalized coverage") +
  theme(legend.position="bottom") +
  sc +
  theme(panel.background = element_rect(fill="gray97")) +
  coord_cartesian(ylim = c(0, 2)) +
  facet_wrap(~ len, ncol=2, labeller = labeller(len = GeneLen_kb.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid"))


## FFT where p=2
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(df[[2]]$ratio_28s_18s),max(df[[2]]$ratio_28s_18s)), name="Ratio 28s/18s")
Ratio28s18s2 <- ggplot(df[[2]], aes(x=normpos, colour=ratio_28s_18s, group=ids)) +
  geom_line(aes(y=df[[2]][, 4]), alpha=1, show.legend = TRUE) +
  geom_hline(yintercept=1, linetype=3) +
  labs(title=sprintf("%s",pairspls[2]),x="Normalized position", y="Mean normalized coverage") +
  theme(legend.position="bottom") +
  sc +
  theme(panel.background = element_rect(fill="gray97")) +
  coord_cartesian(ylim = c(0, 2)) +
  facet_wrap(~ len, ncol=2, labeller = labeller(len = GeneLen_kb.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid"))

sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(df[[2]]$rinvalue),max(df[[2]]$rinvalue)), name="RIN")
RIN2 <- ggplot(df[[2]], aes(x=normpos, colour=rinvalue, group=ids)) +
  geom_line(aes(y=df[[2]][, 4]), alpha=1, show.legend = TRUE) +
  geom_hline(yintercept=1, linetype=3) +
  labs(title=sprintf("%s",pairspls[2]),x="Normalized position", y="Mean normalized coverage") +
  theme(legend.position="bottom") +
  sc +
  theme(panel.background = element_rect(fill="gray97")) +
  coord_cartesian(ylim = c(0, 2)) +
  facet_wrap(~ len, ncol=2, labeller = labeller(len = GeneLen_kb.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid"))

sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(df[[2]]$RatioIntron),max(df[[2]]$RatioIntron)), name="Ratio INTRONIC BASES/CODING BASES")
RatioIntron2 <- ggplot(df[[2]], aes(x=normpos, colour=RatioIntron, group=ids)) +
  geom_line(aes(y=df[[2]][, 4]), alpha=1, show.legend = TRUE) +
  geom_hline(yintercept=1, linetype=3) +
  labs(title=sprintf("%s",pairspls[2]),x="Normalized position", y="Mean normalized coverage") +
  theme(legend.position="bottom") +
  sc +
  theme(panel.background = element_rect(fill="gray97")) +
  coord_cartesian(ylim = c(0, 2)) +
  facet_wrap(~ len, ncol=2, labeller = labeller(len = GeneLen_kb.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid"))


## PET where p=3
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(df[[3]]$ratio_28s_18s),max(df[[3]]$ratio_28s_18s)), name="Ratio 28s/18s")
Ratio28s18s3 <- ggplot(df[[3]], aes(x=normpos, colour=ratio_28s_18s, group=ids)) +
  geom_line(aes(y=df[[3]][, 4]), alpha=1, show.legend = TRUE) +
  geom_hline(yintercept=1, linetype=3) +
  labs(title=sprintf("%s",pairspls[3]),x="Normalized position", y="Mean normalized coverage") +
  theme(legend.position="bottom") +
  sc +
  theme(panel.background = element_rect(fill="gray97")) +
  coord_cartesian(ylim = c(0, 2)) +
  facet_wrap(~ len, ncol=2, labeller = labeller(len = GeneLen_kb.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid"))

sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(df[[3]]$rinvalue),max(df[[3]]$rinvalue)), name="RIN")
RIN3 <- ggplot(df[[3]], aes(x=normpos, colour=rinvalue, group=ids)) +
  geom_line(aes(y=df[[3]][, 4]), alpha=1, show.legend = TRUE) +
  geom_hline(yintercept=1, linetype=3) +
  labs(title=sprintf("%s",pairspls[3]),x="Normalized position", y="Mean normalized coverage") +
  theme(legend.position="bottom") +
  sc +
  theme(panel.background = element_rect(fill="gray97")) +
  coord_cartesian(ylim = c(0, 2)) +
  facet_wrap(~ len, ncol=2, labeller = labeller(len = GeneLen_kb.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid"))

sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(df[[3]]$RatioIntron),max(df[[3]]$RatioIntron)), name="Ratio INTRONIC BASES/CODING BASES")
RatioIntron3 <- ggplot(df[[3]], aes(x=normpos, colour=RatioIntron, group=ids)) +
  geom_line(aes(y=df[[3]][, 4]), alpha=1, show.legend = TRUE) +
  geom_hline(yintercept=1, linetype=3) +
  labs(title=sprintf("%s",pairspls[3]),x="Normalized position", y="Mean normalized coverage") +
  theme(legend.position="bottom") +
  sc +
  theme(panel.background = element_rect(fill="gray97")) +
  coord_cartesian(ylim = c(0, 2)) +
  facet_wrap(~ len, ncol=2, labeller = labeller(len = GeneLen_kb.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid"))


pdf(file="05_PicardMostExpressed_output/PicardMostExpressed_NormTransCoverageByDegra_94ids2.pdf", width=12, height=9,
    bg="white", colormodel="cmyk", paper="letter", onefile=TRUE)
# FFM
grid.arrange(grobs=list(Ratio28s18s1))
grid.arrange(grobs=list(RIN1))
grid.arrange(grobs=list(RatioIntron1))

# FFT
grid.arrange(grobs=list(Ratio28s18s2))
grid.arrange(grobs=list(RIN2))
grid.arrange(grobs=list(RatioIntron2))

# PET
grid.arrange(grobs=list(Ratio28s18s3))
grid.arrange(grobs=list(RIN3))
grid.arrange(grobs=list(RatioIntron3))
dev.off()
