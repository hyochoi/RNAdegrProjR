## --------------------------------------------
##    code for FFT vs. PET mean scatter plot
## --------------------------------------------

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


## ---- 1. Convert FCmat matrix to FCdf dataframe

load("data/FCmat_TMM_SepByPair_log2.rda")
FCmat_TMM_SepByPair1_log2 = FCmat_TMM_SepByPair_log2[,1:26]
FCmat_TMM_SepByPair2_log2 = FCmat_TMM_SepByPair_log2[,27:60]
FCmat_TMM_SepByPair3_log2 = FCmat_TMM_SepByPair_log2[,61:94]

convert_FCmat2FCdf = function(FCmat11) {
  FCmat12 <- t(FCmat11)
  sample <- substr(rownames(FCmat12), start=1, stop=12)
  pair <- substr(rownames(FCmat12), start=14, stop=16)
  FCmat13 <- cbind(sample, pair, FCmat12)

  # Vectorize the counts
  col12 <- repmat(FCmat13[,1:2], (ncol(FCmat13)-2), 1) # (The number of sample ids*The number of genes) x 1
  genes <- matrix(colnames(FCmat13[,3:ncol(FCmat13)]),(ncol(FCmat13)-2),1)
  col30 <- repmat(genes, nrow(FCmat13), 1) # (The number of genes*The number of sample ids) x 1
  col3 <- col30[order(col30),]
  col4 <- as.vector(FCmat13[,3:ncol(FCmat13)]) # (The number of sample ids*The number of genes) x 1
  df <- data.frame(col12, col3, col4)

  # Convert dataframe column from character to numeric
  FCdf = convert_Chr2Numcol(df, 4)
  colnames(FCdf) <- c("sample","pair","geneSymbol","counts")

  return(FCdf)
}

FCdf_TMM_SepByPair1_log2 = convert_FCmat2FCdf(FCmat_TMM_SepByPair1_log2)
FCdf_TMM_SepByPair2_log2 = convert_FCmat2FCdf(FCmat_TMM_SepByPair2_log2)
FCdf_TMM_SepByPair3_log2 = convert_FCmat2FCdf(FCmat_TMM_SepByPair3_log2)
FCdf_TMM_SepByPair_log2 = rbind(FCdf_TMM_SepByPair1_log2, FCdf_TMM_SepByPair2_log2, FCdf_TMM_SepByPair3_log2)


## ---- 2. FFT vs. PET mean scatter plot

FFT1 <- FCdf_TMM_SepByPair_log2[FCdf_TMM_SepByPair_log2$pair=="FFT",]
FFT2 <- aggregate(FFT1$counts, by=list(FFT1$geneSymbol), FUN=mean)
colnames(FFT2) <- c("geneSymbol","FFTmean")

PET1 <- FCdf_TMM_SepByPair_log2[FCdf_TMM_SepByPair_log2$pair=="PET",]
PET2 <- aggregate(PET1$counts, by=list(PET1$geneSymbol), FUN=mean)
colnames(PET2) <- c("geneSymbol","PETmean")

df <- merge(FFT2, PET2, by="geneSymbol", all.x=TRUE, sort=TRUE)

# Add legends
load("data/GeneInfo_KeepTrue2.rda")

ginfo1 <- GeneInfo_KeepTrue2[,c("geneSymbol","median","merged","gene_category","exon.wtpct_gc","subcategory")]
df2 <- merge(df, ginfo1, by="geneSymbol", all.x=TRUE, sort=TRUE)
dim(df2) # 24798     8

# Convert continuous to discrete variables
df3 <- df2
df3$median_grp[df2$median>=0 & df2$median<1000] <- "1_0-1kb"
df3$median_grp[df2$median>=1000 & df2$median<5000] <- "2_1-5kb"
df3$median_grp[df2$median>=5000 & df2$median<10000] <- "3_5-10kb"
df3$median_grp[df2$median>=10000] <- "4_10+kb"

df3$merged_grp[df3$merged>=0 & df3$merged<1000] <- "1_0-1kb"
df3$merged_grp[df3$merged>=1000 & df3$merged<5000] <- "2_1-5kb"
df3$merged_grp[df3$merged>=5000 & df3$merged<10000] <- "3_5-10kb"
df3$merged_grp[df3$merged>=10000] <- "4_10+kb"

quan <- quantile(df3$exon.wtpct_gc, c(0.25, 0.50, 0.75))
df3$exon.wtpct_gc_grp[df3$exon.wtpct_gc<quan[1]] <- "Min-Q1"
df3$exon.wtpct_gc_grp[df3$exon.wtpct_gc>=quan[1] & df3$exon.wtpct_gc<quan[2]] <- "Q1-Q2"
df3$exon.wtpct_gc_grp[df3$exon.wtpct_gc>=quan[2] & df3$exon.wtpct_gc<quan[3]] <- "Q2-Q3"
df3$exon.wtpct_gc_grp[df3$exon.wtpct_gc>=quan[3]] <- "Q3-Max"

# Find ACTB, GAPDH, MALAT1, NEAT1
ACTB <- as.matrix(df3[df3$geneSymbol=="ACTB",c("FFTmean","PETmean")])
GAPDH <- as.matrix(df3[df3$geneSymbol=="GAPDH",c("FFTmean","PETmean")])
MALAT1 <- as.matrix(df3[df3$geneSymbol=="MALAT1",c("FFTmean","PETmean")])
NEAT1 <- as.matrix(df3[df3$geneSymbol=="NEAT1",c("FFTmean","PETmean")])


## Legend: subcategory
# Remove NA, TEC (to be experimentally confirmed) in the subcategory figures
s1 <- ggplot(df3[!(is.na(df3$subcategory)) & df3$subcategory!="TEC",], aes(x=FFTmean, y=PETmean, color=subcategory)) +
  geom_point(alpha=0.2) + geom_rug(alpha=0.1) +
  scale_colour_brewer(palette="Dark2") +
  labs(title="FFT vs. PET mean counts of log2(TMM+1) by 24,798 genes after filtering; \nremoved NA, TEC (to be experimentally confirmed)",x="FFT mean", y="PET mean") +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(title="Subcategory",nrow=1,byrow=TRUE)) +
  geom_abline(intercept=0, slope=1, linetype=3, color="black", size=1.2) +
  geom_smooth(method=lm, se=FALSE) +
  theme(panel.background = element_rect(fill="gray97")) +
  geom_point(aes(x=ACTB[1], y=ACTB[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=GAPDH[1], y=GAPDH[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=MALAT1[1], y=MALAT1[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=NEAT1[1], y=NEAT1[2]), colour="black", shape=16, size=2) +
  geom_text(aes(label=ifelse(geneSymbol=="ACTB"|geneSymbol=="GAPDH"|geneSymbol=="MALAT1"|geneSymbol=="NEAT1",as.character(geneSymbol),'')), hjust=-0.1, vjust=1)

s2 <- ggplot(df3[!(is.na(df3$subcategory)) & df3$subcategory!="TEC",], aes(x=FFTmean, y=PETmean, color=subcategory)) +
  geom_point(alpha=0.2) + geom_rug(alpha=0.1) +
  scale_colour_brewer(palette="Dark2") +
  labs(title="FFT vs. PET mean counts of log2(TMM+1) by 24,798 genes after filtering; \nremoved NA, TEC (to be experimentally confirmed)",x="FFT mean", y="PET mean") +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(title="Subcategory",nrow=1,byrow=TRUE)) +
  geom_abline(intercept=0, slope=1, linetype=3, color="black", size=1.2) +
  geom_smooth(method=lm, se=FALSE) +
  theme(panel.background = element_rect(fill="gray97")) +
  facet_wrap(~ subcategory, ncol=2) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid")) +
  geom_point(aes(x=ACTB[1], y=ACTB[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=GAPDH[1], y=GAPDH[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=MALAT1[1], y=MALAT1[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=NEAT1[1], y=NEAT1[2]), colour="black", shape=16, size=2) +
  geom_text(aes(label=ifelse(geneSymbol=="ACTB"|geneSymbol=="GAPDH"|geneSymbol=="MALAT1"|geneSymbol=="NEAT1",as.character(geneSymbol),'')), hjust=-0.1, vjust=1)


## Legend: median
# Continuous legends
myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))

sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(log10(df3$median)),max(log10(df3$median))), name="log10(Median)")
m1 <- ggplot(df3, aes(x=FFTmean, y=PETmean, colour=log10(median))) +
  geom_point(alpha=0.4) + geom_rug(alpha=0.1) +
  labs(title="FFT vs. PET mean counts of log2(TMM+1) by 24,798 genes after filtering",x="FFT mean", y="PET mean") +
  theme(legend.position="bottom") +
  sc +
  geom_abline(intercept=0, slope=1, linetype=3, color="black", size=1.2) +
  geom_smooth(method=lm, se=FALSE) +
  theme(panel.background = element_rect(fill="gray97")) +
  geom_point(aes(x=ACTB[1], y=ACTB[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=GAPDH[1], y=GAPDH[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=MALAT1[1], y=MALAT1[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=NEAT1[1], y=NEAT1[2]), colour="black", shape=16, size=2) +
  geom_text(aes(label=ifelse(geneSymbol=="ACTB"|geneSymbol=="GAPDH"|geneSymbol=="MALAT1"|geneSymbol=="NEAT1",as.character(geneSymbol),'')), hjust=-0.1, vjust=1)

# New facet label names
median_grp.labs <- c("0-1kb", "1-5kb", "5-10kb", "10+kb")
names(median_grp.labs) <- c("1_0-1kb", "2_1-5kb", "3_5-10kb", "4_10+kb")

m2 <- ggplot(df3, aes(x=FFTmean, y=PETmean, colour=median_grp)) +
  geom_point(alpha=0.2) + geom_rug(alpha=0.1) +
  scale_colour_brewer(palette="Dark2", labels=c("0-1kb", "1-5kb", "5-10kb", "10+kb")) +
  labs(title="FFT vs. PET mean counts of log2(TMM+1) by 24,798 genes after filtering",x="FFT mean", y="PET mean") +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(title="Median",nrow=1,byrow=TRUE)) +
  geom_abline(intercept=0, slope=1, linetype=3, color="black", size=1.2) +
  geom_smooth(method=lm, se=FALSE) +
  theme(panel.background = element_rect(fill="gray97")) +
  facet_wrap(~ median_grp, ncol=2, labeller = labeller(median_grp = median_grp.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid")) +
  geom_point(aes(x=ACTB[1], y=ACTB[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=GAPDH[1], y=GAPDH[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=MALAT1[1], y=MALAT1[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=NEAT1[1], y=NEAT1[2]), colour="black", shape=16, size=2) +
  geom_text(aes(label=ifelse(geneSymbol=="ACTB"|geneSymbol=="GAPDH"|geneSymbol=="MALAT1"|geneSymbol=="NEAT1",as.character(geneSymbol),'')), hjust=-0.1, vjust=1)


## Legend: merged
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(log10(df3$merged)),max(log10(df3$merged))), name="log10(Merged)")
r1 <- ggplot(df3, aes(x=FFTmean, y=PETmean, colour=log10(merged))) +
  geom_point(alpha=0.4) + geom_rug(alpha=0.1) +
  labs(title="FFT vs. PET mean counts of log2(TMM+1) by 24,798 genes after filtering",x="FFT mean", y="PET mean") +
  theme(legend.position="bottom") +
  sc +
  geom_abline(intercept=0, slope=1, linetype=3, color="black", size=1.2) +
  geom_smooth(method=lm, se=FALSE) +
  theme(panel.background = element_rect(fill="gray97")) +
  geom_point(aes(x=ACTB[1], y=ACTB[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=GAPDH[1], y=GAPDH[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=MALAT1[1], y=MALAT1[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=NEAT1[1], y=NEAT1[2]), colour="black", shape=16, size=2) +
  geom_text(aes(label=ifelse(geneSymbol=="ACTB"|geneSymbol=="GAPDH"|geneSymbol=="MALAT1"|geneSymbol=="NEAT1",as.character(geneSymbol),'')), hjust=-0.1, vjust=1)

# New facet label names
merged_grp.labs <- c("0-1kb", "1-5kb", "5-10kb", "10+kb")
names(merged_grp.labs) <- c("1_0-1kb", "2_1-5kb", "3_5-10kb", "4_10+kb")

r2 <- ggplot(df3, aes(x=FFTmean, y=PETmean, colour=merged_grp)) +
  geom_point(alpha=0.2) + geom_rug(alpha=0.1) +
  scale_colour_brewer(palette="Dark2", labels=c("0-1kb", "1-5kb", "5-10kb", "10+kb")) +
  labs(title="FFT vs. PET mean counts of log2(TMM+1) by 24,798 genes after filtering",x="FFT mean", y="PET mean") +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(title="Merged",nrow=1,byrow=TRUE)) +
  geom_abline(intercept=0, slope=1, linetype=3, color="black", size=1.2) +
  geom_smooth(method=lm, se=FALSE) +
  theme(panel.background = element_rect(fill="gray97")) +
  facet_wrap(~ merged_grp, ncol=2, labeller = labeller(merged_grp = merged_grp.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid")) +
  geom_point(aes(x=ACTB[1], y=ACTB[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=GAPDH[1], y=GAPDH[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=MALAT1[1], y=MALAT1[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=NEAT1[1], y=NEAT1[2]), colour="black", shape=16, size=2) +
  geom_text(aes(label=ifelse(geneSymbol=="ACTB"|geneSymbol=="GAPDH"|geneSymbol=="MALAT1"|geneSymbol=="NEAT1",as.character(geneSymbol),'')), hjust=-0.1, vjust=1)

r3 <- ggplot(df3[!(is.na(df3$subcategory)) & df3$subcategory!="TEC",], aes(x=FFTmean, y=PETmean, colour=merged_grp)) +
  geom_point(alpha=0.2) + geom_rug(alpha=0.1) +
  scale_colour_brewer(palette="Dark2", labels=c("0-1kb", "1-5kb", "5-10kb", "10+kb")) +
  labs(title="FFT vs. PET mean counts of log2(TMM+1) by 24,798 genes after filtering; \nremoved NA, TEC (to be experimentally confirmed)",x="FFT mean", y="PET mean") +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(title="Merged",nrow=1,byrow=TRUE)) +
  geom_abline(intercept=0, slope=1, linetype=3, color="black", size=1.2) +
  geom_smooth(method=lm, se=FALSE) +
  theme(panel.background = element_rect(fill="gray97")) +
  facet_grid(subcategory ~ merged_grp, labeller = labeller(merged_grp = merged_grp.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid")) +
  geom_point(aes(x=ACTB[1], y=ACTB[2]), colour="black", shape=16, size=1.5) +
  geom_point(aes(x=GAPDH[1], y=GAPDH[2]), colour="black", shape=16, size=1.5) +
  geom_point(aes(x=MALAT1[1], y=MALAT1[2]), colour="black", shape=16, size=1.5) +
  geom_point(aes(x=NEAT1[1], y=NEAT1[2]), colour="black", shape=16, size=1.5) +
  geom_text(size=2.5, aes(label=ifelse(geneSymbol=="ACTB"|geneSymbol=="GAPDH"|geneSymbol=="MALAT1"|geneSymbol=="NEAT1",as.character(geneSymbol),'')), hjust=-0.1, vjust=1)


## Legend: exon.wtpct_gc
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(df3$exon.wtpct_gc),max(df3$exon.wtpct_gc)), name="Exon weighted percentage of GC")
e1 <- ggplot(df3, aes(x=FFTmean, y=PETmean, colour=exon.wtpct_gc)) +
  geom_point(alpha=0.4) + geom_rug(alpha=0.1) +
  labs(title="FFT vs. PET mean counts of log2(TMM+1) by 24,798 genes after filtering",x="FFT mean", y="PET mean") +
  theme(legend.position="bottom") +
  sc +
  geom_abline(intercept=0, slope=1, linetype=3, color="black", size=1.2) +
  geom_smooth(method=lm, se=FALSE) +
  theme(panel.background = element_rect(fill="gray97")) +
  geom_point(aes(x=ACTB[1], y=ACTB[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=GAPDH[1], y=GAPDH[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=MALAT1[1], y=MALAT1[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=NEAT1[1], y=NEAT1[2]), colour="black", shape=16, size=2) +
  geom_text(aes(label=ifelse(geneSymbol=="ACTB"|geneSymbol=="GAPDH"|geneSymbol=="MALAT1"|geneSymbol=="NEAT1",as.character(geneSymbol),'')), hjust=-0.1, vjust=1)

e2 <- ggplot(df3, aes(x=FFTmean, y=PETmean, colour=exon.wtpct_gc_grp)) +
  geom_point(alpha=0.2) + geom_rug(alpha=0.1) +
  scale_colour_brewer(palette="Dark2") +
  labs(title="FFT vs. PET mean counts of log2(TMM+1) by 24,798 genes after filtering",x="FFT mean", y="PET mean") +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(title="Exon weighted percentage of GC",nrow=1,byrow=TRUE)) +
  geom_abline(intercept=0, slope=1, linetype=3, color="black", size=1.2) +
  geom_smooth(method=lm, se=FALSE) +
  theme(panel.background = element_rect(fill="gray97")) +
  facet_wrap(~ exon.wtpct_gc_grp, ncol=2) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid")) +
  geom_point(aes(x=ACTB[1], y=ACTB[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=GAPDH[1], y=GAPDH[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=MALAT1[1], y=MALAT1[2]), colour="black", shape=16, size=2) +
  geom_point(aes(x=NEAT1[1], y=NEAT1[2]), colour="black", shape=16, size=2) +
  geom_text(aes(label=ifelse(geneSymbol=="ACTB"|geneSymbol=="GAPDH"|geneSymbol=="MALAT1"|geneSymbol=="NEAT1",as.character(geneSymbol),'')), hjust=-0.1, vjust=1)

e3 <- ggplot(df3[!(is.na(df3$subcategory)) & df3$subcategory!="TEC",], aes(x=FFTmean, y=PETmean, colour=exon.wtpct_gc_grp)) +
  geom_point(alpha=0.2) + geom_rug(alpha=0.1) +
  scale_colour_brewer(palette="Dark2") +
  labs(title="FFT vs. PET mean counts of log2(TMM+1) by 24,798 genes after filtering; \nremoved NA, TEC (to be experimentally confirmed)",x="FFT mean", y="PET mean") +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(title="Exon weighted percentage of GC",nrow=1,byrow=TRUE)) +
  geom_abline(intercept=0, slope=1, linetype=3, color="black", size=1.2) +
  geom_smooth(method=lm, se=FALSE) +
  theme(panel.background = element_rect(fill="gray97")) +
  facet_grid(subcategory ~ exon.wtpct_gc_grp) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid")) +
  geom_point(aes(x=ACTB[1], y=ACTB[2]), colour="black", shape=16, size=1.5) +
  geom_point(aes(x=GAPDH[1], y=GAPDH[2]), colour="black", shape=16, size=1.5) +
  geom_point(aes(x=MALAT1[1], y=MALAT1[2]), colour="black", shape=16, size=1.5) +
  geom_point(aes(x=NEAT1[1], y=NEAT1[2]), colour="black", shape=16, size=1.5) +
  geom_text(size=2.5, aes(label=ifelse(geneSymbol=="ACTB"|geneSymbol=="GAPDH"|geneSymbol=="MALAT1"|geneSymbol=="NEAT1",as.character(geneSymbol),'')), hjust=-0.1, vjust=1)

e4 <- ggplot(df3, aes(x=FFTmean, y=PETmean, colour=exon.wtpct_gc_grp)) +
  geom_point(alpha=0.2) + geom_rug(alpha=0.1) +
  scale_colour_brewer(palette="Dark2") +
  labs(title="FFT vs. PET mean counts of log2(TMM+1) by 24,798 genes after filtering",x="FFT mean", y="PET mean") +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(title="Exon weighted percentage of GC",nrow=1,byrow=TRUE)) +
  geom_abline(intercept=0, slope=1, linetype=3, color="black", size=1.2) +
  geom_smooth(method=lm, se=FALSE) +
  theme(panel.background = element_rect(fill="gray97")) +
  facet_grid(merged_grp ~ exon.wtpct_gc_grp, labeller = labeller(merged_grp = merged_grp.labs)) +
  theme(
    strip.text.x = element_text(size = 12, color = "black"),
    strip.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(color="NA", fill="white", linewidth=1, linetype="solid")) +
  geom_point(aes(x=ACTB[1], y=ACTB[2]), colour="black", shape=16, size=1.5) +
  geom_point(aes(x=GAPDH[1], y=GAPDH[2]), colour="black", shape=16, size=1.5) +
  geom_point(aes(x=MALAT1[1], y=MALAT1[2]), colour="black", shape=16, size=1.5) +
  geom_point(aes(x=NEAT1[1], y=NEAT1[2]), colour="black", shape=16, size=1.5) +
  geom_text(size=2.5, aes(label=ifelse(geneSymbol=="ACTB"|geneSymbol=="GAPDH"|geneSymbol=="MALAT1"|geneSymbol=="NEAT1",as.character(geneSymbol),'')), hjust=-0.1, vjust=1)


pdf(file="01_DataPreprocessing_output/FragmentCounts_norm_94ids_AllCategory_filterByExpr/FFTPETmeanScatterplot_TMMlog2_filterByExpr4.pdf", width=12, height=9,
    bg="white", colormodel="cmyk", paper="letter", onefile=TRUE)
grid.arrange(grobs=list(s1))
grid.arrange(grobs=list(s2))
grid.arrange(grobs=list(m1))
grid.arrange(grobs=list(m2))
grid.arrange(grobs=list(r1))
grid.arrange(grobs=list(r2))
grid.arrange(grobs=list(r3))
grid.arrange(grobs=list(e1))
grid.arrange(grobs=list(e2))
grid.arrange(grobs=list(e3))
grid.arrange(grobs=list(e4))
dev.off()
