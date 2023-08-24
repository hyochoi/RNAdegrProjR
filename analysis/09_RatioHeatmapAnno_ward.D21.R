## ---------------------------------------------------
##    code for ratio heatmap with ward.D clustering
## ---------------------------------------------------

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
  library(ComplexHeatmap) # For heatmap
  library(circlize) # For color interpolation
  library(colortools) # For color wheel
})


## ---- 1. PET - FFT

load("data/FCmat_TMM_SepByPair_log2.rda")

# Separate columns by pair
pair <- substr(colnames(FCmat_TMM_SepByPair_log2), start=14, stop=16)
keepcol <- c("FFM","FFT","PET")

FFTcol <- grep(pattern=keepcol[2],x=pair)
PETcol <- grep(pattern=keepcol[3],x=pair)

FFTmat <- FCmat_TMM_SepByPair_log2[,c(FFTcol)]
PETmat <- FCmat_TMM_SepByPair_log2[,c(PETcol)]

# Sort by colnames (patients)
FFTmat1 <- FFTmat[, order(colnames(FFTmat))]
PETmat1 <- PETmat[, order(colnames(PETmat))]

# Sort by rownames (genes)
FFTmat2 <- FFTmat1[order(rownames(FFTmat1)), ]
PETmat2 <- PETmat1[order(rownames(PETmat1)), ]

# PET - FFT
Ratiomat = PETmat2-FFTmat2
colnames(Ratiomat) <- substr(colnames(PETmat2), start=1, stop=12)
dim(Ratiomat) # 24798    34


## ---- 2. Heatmap

## ---- 2.1 Gene info annotation

load("data-raw/FCdf_TMM_SepByPair_log2.RData") # from 01_FFTPETmeanScatterplot_FCViolinplot.R

FFT1 <- FCdf_TMM_SepByPair_log2[FCdf_TMM_SepByPair_log2$pair=="FFT",]
FFT2 <- aggregate(FFT1$counts, by=list(FFT1$geneSymbol), FUN=mean)
colnames(FFT2) <- c("geneSymbol","FFTmean")

PET1 <- FCdf_TMM_SepByPair_log2[FCdf_TMM_SepByPair_log2$pair=="PET",]
PET2 <- aggregate(PET1$counts, by=list(PET1$geneSymbol), FUN=mean)
colnames(PET2) <- c("geneSymbol","PETmean")

df <- merge(FFT2, PET2, by="geneSymbol", all.x=TRUE, sort=TRUE)

# RatioPETFFT
df$RatioPETFFT = df$PETmean-df$FFTmean

# Sort by RatioPETFFT
df3 <- df[order(df$RatioPETFFT, decreasing=FALSE), ]
x <- c(1:nrow(df3))
df4 <- cbind(x, df3)

# Add gene info
load("data/GeneInfo_KeepTrue2.rda")

df6 <- merge(df4, GeneInfo_KeepTrue2, by="geneSymbol", all.x=TRUE, sort=TRUE)
rownames(df6) <- df6$geneSymbol

# Sort by RatioPETFFT
df7 <- df6[order(df6$RatioPETFFT, decreasing=FALSE), ]
df8 <- df7[,c("RatioPETFFT","FFTmean","PETmean","merged","exon.wtpct_gc","intron.wtpct_gc","3UTR.meanlen","5UTR.meanlen","subcategory")]

# Convert a continuous variable to eCDF values
mt <- t(data.matrix(df8))

FFTmeandf = convert_Cont2eCDF.gene(mt[c("FFTmean"),],"FFTmean")
PETmeandf = convert_Cont2eCDF.gene(mt[c("PETmean"),],"PETmean")
mergeddf = convert_Cont2eCDF.gene(mt[c("merged"),],"merged")
meanlen.3UTRdf = convert_Cont2eCDF.gene(mt[c("3UTR.meanlen"),],"3UTR.meanlen")
meanlen.5UTRdf = convert_Cont2eCDF.gene(mt[c("5UTR.meanlen"),],"5UTR.meanlen")
exon.wtpct_gcdf = convert_Cont2eCDF.gene(mt[c("exon.wtpct_gc"),],"exon.wtpct_gc")
intron.wtpct_gcdf = convert_Cont2eCDF.gene(mt[c("intron.wtpct_gc"),],"intron.wtpct_gc")

geneSymbol <- as.matrix(rownames(df8))
df9 <- cbind(geneSymbol, df8)
df10 <- merge(df9, FFTmeandf, by="geneSymbol", all.x=TRUE)
df11 <- merge(df10, PETmeandf, by="geneSymbol", all.x=TRUE)
df12 <- merge(df11, mergeddf, by="geneSymbol", all.x=TRUE)
df13 <- merge(df12, meanlen.3UTRdf, by="geneSymbol", all.x=TRUE)
df14 <- merge(df13, meanlen.5UTRdf, by="geneSymbol", all.x=TRUE)
df15 <- merge(df14, exon.wtpct_gcdf, by="geneSymbol", all.x=TRUE)
df16 <- merge(df15, intron.wtpct_gcdf, by="geneSymbol", all.x=TRUE)

# Sort by RatioPETFFT
df17 <- df16[order(df16$RatioPETFFT, decreasing=FALSE), ]
rownames(df17) <- c(df17$geneSymbol)

# Sort by gene names
GeneAnnomat <- df17[order(rownames(df17)), ]
dim(GeneAnnomat) # 24798    17


## rowAnnotation
# Find complementary colors
colwheel <- wheel("#E41A1C", num = 32)
[1] "#E41A1C" "#E43E1A" "#E4641A" "#E48A1A" "#E4B01A" "#E4D51A" "#CDE41A" "#A7E41A" "#81E41A" "#5BE41A" "#35E41A" "#1AE425" "#1AE44B" "#1AE470" "#1AE496" "#1AE4BC"
    "#1AE4E2" "#1AC0E4" "#1A9AE4" "#1A74E4" "#1A4EE4" "#1A29E4" "#311AE4" "#571AE4" "#7D1AE4" "#A31AE4" "#C91AE4" "#E41AD9" "#E41AB3" "#E41A8E" "#E41A68" "#E41A42"

col_sc = c("lncRNA" = "#E41A1C", "non-coding RNA" = "#377EB8", "protein_coding" = "#4DAF4A", "pseudogene" = "#984EA3", "TEC" = "#FF7F00")
lgd1 = Legend(labels = c("LncRNA","Non-coding RNA","Protein coding","Pseudogene","TEC"), title = "Subcategory", legend_gp = gpar(fill = col_sc), title_position = "lefttop", labels_gp = gpar(fontsize = 8))

col_fun_prop2 = colorRamp2(c(0, 0.5, 1), c(colwheel[17], "white", colwheel[1]))
lgd2 = Legend(col_fun = col_fun_prop2, title = "FFT mean", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
              labels = c(round(min(GeneAnnomat$FFTmean),1), round(quantile(GeneAnnomat$FFTmean, 0.5),1), round(max(GeneAnnomat$FFTmean),1)))

col_fun_prop3 = colorRamp2(c(0, 0.5, 1), c(colwheel[18], "white", colwheel[2]))
lgd3 = Legend(col_fun = col_fun_prop3, title = "PET mean", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
              labels = c(round(min(GeneAnnomat$PETmean),1), round(quantile(GeneAnnomat$PETmean, 0.5),1), round(max(GeneAnnomat$PETmean),1)))

col_fun_prop4 = colorRamp2(c(0, 0.5, 1), c(colwheel[19], "white", colwheel[3]))
lgd4 = Legend(col_fun = col_fun_prop4, title = "Gene length", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
              labels = c(round(min(GeneAnnomat$merged),1), round(quantile(GeneAnnomat$merged, 0.5),1), round(max(GeneAnnomat$merged),1)))

col_fun_prop5 = colorRamp2(c(0, 0.5, 1), c(colwheel[20], "white", colwheel[4]))
lgd5 = Legend(col_fun = col_fun_prop5, title = "3UTR mean length", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
              labels = c(round(min(GeneAnnomat[,c("3UTR.meanlen")],na.rm=TRUE),1), round(quantile(GeneAnnomat[,c("3UTR.meanlen")], 0.5,na.rm=TRUE),1), round(max(GeneAnnomat[,c("3UTR.meanlen")],na.rm=TRUE),1)))

col_fun_prop6 = colorRamp2(c(0, 0.5, 1), c(colwheel[21], "white", colwheel[5]))
lgd6 = Legend(col_fun = col_fun_prop6, title = "5UTR mean length", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
              labels = c(round(min(GeneAnnomat[,c("5UTR.meanlen")],na.rm=TRUE),1), round(quantile(GeneAnnomat[,c("5UTR.meanlen")], 0.5,na.rm=TRUE),1), round(max(GeneAnnomat[,c("5UTR.meanlen")],na.rm=TRUE),1)))

col_fun_prop7 = colorRamp2(c(0, 0.5, 1), c(colwheel[22], "white", colwheel[6]))
lgd7 = Legend(col_fun = col_fun_prop7, title = "Exon GC percentage", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
              labels = c(round(min(GeneAnnomat$exon.wtpct_gc),1), round(quantile(GeneAnnomat$exon.wtpct_gc, 0.5),1), round(max(GeneAnnomat$exon.wtpct_gc),1)))

col_fun_prop8 = colorRamp2(c(0, 0.5, 1), c(colwheel[23], "white", colwheel[7]))
lgd8 = Legend(col_fun = col_fun_prop8, title = "Intron GC percentage", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
              labels = c(round(min(GeneAnnomat$intron.wtpct_gc,na.rm=TRUE),1), round(quantile(GeneAnnomat$intron.wtpct_gc, 0.5,na.rm=TRUE),1), round(max(GeneAnnomat$intron.wtpct_gc,na.rm=TRUE),1)))

# Find ACTB, GAPDH, MALAT1, NEAT1 (exact match)
ACTB <- grep(pattern="^ACTB$",x=rownames(Ratiomat)) # 3575
GAPDH <- grep(pattern="^GAPDH$",x=rownames(Ratiomat)) # 10780
MALAT1 <- grep(pattern="^MALAT1$",x=rownames(Ratiomat)) # 14243
NEAT1 <- grep(pattern="^NEAT1$",x=rownames(Ratiomat)) # 15713

row_ha5 = rowAnnotation(mark = anno_mark(at = c(ACTB, GAPDH, MALAT1, NEAT1), labels = c("ACTB", "GAPDH", "MALAT1", "NEAT1"), side = "left",
                                         labels_gp = gpar(fontsize = 8, col = c(rep("#4DAF4A", 2), rep("#E41A1C", 2))), labels_rot = 90),
                        Subcategory = GeneAnnomat$subcategory,
                        "PET mean − FFT mean" = anno_points(GeneAnnomat$RatioPETFFT, size = unit(0.5, "mm"),
                                                            gp = gpar(col = ifelse(GeneAnnomat$RatioPETFFT > 0, "red", "blue")),
                                                            axis_param = list(labels_rot = 0)),
                        "FFT mean" = GeneAnnomat$FFTmean_ecdf,
                        "PET mean" = GeneAnnomat$PETmean_ecdf,
                        "Gene length" = GeneAnnomat$merged_ecdf,
                        "3UTR mean length" = GeneAnnomat[,c("3UTR.meanlen_ecdf")],
                        "5UTR mean length" = GeneAnnomat[,c("5UTR.meanlen_ecdf")],
                        "Exon GC percentage" = GeneAnnomat$exon.wtpct_gc_ecdf,
                        "Intron GC percentage" = GeneAnnomat$intron.wtpct_gc_ecdf,
                        annotation_name_gp= gpar(fontsize = 8), annotation_name_rot = c(90),
                        col = list(Subcategory = col_sc,
                                   "FFT mean" = col_fun_prop2,
                                   "PET mean" = col_fun_prop3,
                                   "Gene length" = col_fun_prop4,
                                   "3UTR mean length" = col_fun_prop5,
                                   "5UTR mean length" = col_fun_prop6,
                                   "Exon GC percentage" = col_fun_prop7,
                                   "Intron GC percentage" = col_fun_prop8), show_legend = FALSE)


## ---- 2.2 Sample info annotation

# Sample information table
load("data/SampleInfo.rda")

# keep 94 ids
removespls <- c("TCGA-A6-2674-FFM",
                "TCGA-A6-2684-FFM",
                "TCGA-A6-3809-FFM",
                "TCGA-A6-3810-FFM",
                "TCGA-BK-A0CA-FFM",
                "TCGA-BK-A0CC-FFM",
                "TCGA-BK-A139-FFM",
                "TCGA-BK-A26L-FFM")
removerownum = matrix(0, 8, 1)
for (r in 1:8){
  removerownum[r,] = grep(pattern=removespls[r],x=SampleInfo$NewSampleId)
}
SampleInfo1 <- SampleInfo[-removerownum,]
dim(SampleInfo1) # 94 39

# Keep PET samples
pair2 <- substr(SampleInfo1$NewSampleId, start=14, stop=16)
SampleInfo2 <- cbind(pair2, SampleInfo1)
SampleInfo5 <- SampleInfo2[SampleInfo2$pair2=="PET",]

# Convert a continuous variable to eCDF values
rownames(SampleInfo5) <- c(SampleInfo5$NewSampleId)
st <- t(data.matrix(SampleInfo5))

PCT_CODING_BASESdf = convert_Cont2eCDF.spl(st[c("PCT_CODING_BASES"),],"PCT_CODING_BASES")
PCT_INTRONIC_BASESdf = convert_Cont2eCDF.spl(st[c("PCT_INTRONIC_BASES"),],"PCT_INTRONIC_BASES")
PCT_INTERGENIC_BASESdf = convert_Cont2eCDF.spl(st[c("PCT_INTERGENIC_BASES"),],"PCT_INTERGENIC_BASES")
PCT_MRNA_BASESdf = convert_Cont2eCDF.spl(st[c("PCT_MRNA_BASES"),],"PCT_MRNA_BASES")
MEDIAN_CV_COVERAGEdf = convert_Cont2eCDF.spl(st[c("MEDIAN_CV_COVERAGE"),],"MEDIAN_CV_COVERAGE")
MEDIAN_3PRIME_BIASdf = convert_Cont2eCDF.spl(st[c("MEDIAN_3PRIME_BIAS"),],"MEDIAN_3PRIME_BIAS")
MEDIAN_5PRIME_BIASdf = convert_Cont2eCDF.spl(st[c("MEDIAN_5PRIME_BIAS"),],"MEDIAN_5PRIME_BIAS")
RatioIntrondf = convert_Cont2eCDF.spl(st[c("RatioIntron"),],"RatioIntron")
rinvaluedf = convert_Cont2eCDF.spl(st[c("rinvalue"),],"rinvalue")

SampleInfo6 <- merge(SampleInfo5, PCT_CODING_BASESdf, by="NewSampleId", all.x=TRUE)
SampleInfo7 <- merge(SampleInfo6, PCT_INTRONIC_BASESdf, by="NewSampleId", all.x=TRUE)
SampleInfo8 <- merge(SampleInfo7, PCT_INTERGENIC_BASESdf, by="NewSampleId", all.x=TRUE)
SampleInfo9 <- merge(SampleInfo8, PCT_MRNA_BASESdf, by="NewSampleId", all.x=TRUE)
SampleInfo10 <- merge(SampleInfo9, MEDIAN_CV_COVERAGEdf, by="NewSampleId", all.x=TRUE)
SampleInfo11 <- merge(SampleInfo10, MEDIAN_3PRIME_BIASdf, by="NewSampleId", all.x=TRUE)
SampleInfo12 <- merge(SampleInfo11, MEDIAN_5PRIME_BIASdf, by="NewSampleId", all.x=TRUE)
SampleInfo13 <- merge(SampleInfo12, RatioIntrondf, by="NewSampleId", all.x=TRUE)
SampleInfo14 <- merge(SampleInfo13, rinvaluedf, by="NewSampleId", all.x=TRUE)

# Sort by patients names
SplAnnomat <- SampleInfo14[order(SampleInfo14$NewSampleId), ]
dim(SplAnnomat) # 34 49


## HeatmapAnnotation
col_ts = c("BLCA" = "#66C2A5", "BRCA" = "#FC8D62", "COAD" = "#8DA0CB", "KIRC" = "#E78AC3", "LUAD" = "#A6D854", "UCEC" = "#FFD92F")
lgd18 = Legend(labels = c("BLCA", "BRCA", "COAD", "KIRC", "LUAD", "UCEC"),
               title = "Tissue", legend_gp = gpar(fill = col_ts), title_position = "lefttop", labels_gp = gpar(fontsize = 8))

col_fun_prop9 = colorRamp2(c(0, 0.5, 1), c(colwheel[24], "white", colwheel[8]))
lgd9 = Legend(col_fun = col_fun_prop9, title = "PET PCT_CODING_BASES", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
              labels = c(round(min(SplAnnomat$PCT_CODING_BASES),1), round(quantile(SplAnnomat$PCT_CODING_BASES, 0.5),1), round(max(SplAnnomat$PCT_CODING_BASES),1)))

col_fun_prop10 = colorRamp2(c(0, 0.5, 1), c(colwheel[25], "white", colwheel[9]))
lgd10 = Legend(col_fun = col_fun_prop10, title = "PET PCT_INTRONIC_BASES", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
               labels = c(round(min(SplAnnomat$PCT_INTRONIC_BASES),1), round(quantile(SplAnnomat$PCT_INTRONIC_BASES, 0.5),1), round(max(SplAnnomat$PCT_INTRONIC_BASES),1)))

col_fun_prop11 = colorRamp2(c(0, 0.5, 1), c(colwheel[26], "white", colwheel[10]))
lgd11 = Legend(col_fun = col_fun_prop11, title = "PET PCT_INTERGENIC_BASES", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
               labels = c(round(min(SplAnnomat$PCT_INTERGENIC_BASES),1), round(quantile(SplAnnomat$PCT_INTERGENIC_BASES, 0.5),1), round(max(SplAnnomat$PCT_INTERGENIC_BASES),1)))

col_fun_prop12 = colorRamp2(c(0, 0.5, 1), c(colwheel[27], "white", colwheel[11]))
lgd12 = Legend(col_fun = col_fun_prop12, title = "PET PCT_MRNA_BASES", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
               labels = c(round(min(SplAnnomat$PCT_MRNA_BASES),1), round(quantile(SplAnnomat$PCT_MRNA_BASES, 0.5),1), round(max(SplAnnomat$PCT_MRNA_BASES),1)))

col_fun_prop13 = colorRamp2(c(0, 0.5, 1), c(colwheel[28], "white", colwheel[12]))
lgd13 = Legend(col_fun = col_fun_prop13, title = "PET MEDIAN_CV_COVERAGE", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
               labels = c(round(min(SplAnnomat$MEDIAN_CV_COVERAGE),1), round(quantile(SplAnnomat$MEDIAN_CV_COVERAGE, 0.5),1), round(max(SplAnnomat$MEDIAN_CV_COVERAGE),1)))

col_fun_prop14 = colorRamp2(c(0, 0.5, 1), c(colwheel[29], "white", colwheel[13]))
lgd14 = Legend(col_fun = col_fun_prop14, title = "PET MEDIAN_3PRIME_BIAS", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
               labels = c(round(min(SplAnnomat$MEDIAN_3PRIME_BIAS),1), round(quantile(SplAnnomat$MEDIAN_3PRIME_BIAS, 0.5),1), round(max(SplAnnomat$MEDIAN_3PRIME_BIAS),1)))

col_fun_prop15 = colorRamp2(c(0, 0.5, 1), c(colwheel[30], "white", colwheel[14]))
lgd15 = Legend(col_fun = col_fun_prop15, title = "PET MEDIAN_5PRIME_BIAS", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
               labels = c(round(min(SplAnnomat$MEDIAN_5PRIME_BIAS),1), round(quantile(SplAnnomat$MEDIAN_5PRIME_BIAS, 0.5),1), round(max(SplAnnomat$MEDIAN_5PRIME_BIAS),1)))

col_fun_prop16 = colorRamp2(c(0, 0.5, 1), c(colwheel[31], "white", colwheel[15]))
lgd16 = Legend(col_fun = col_fun_prop16, title = "PET RatioIntron", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
               labels = c(round(min(SplAnnomat$RatioIntron),1), round(quantile(SplAnnomat$RatioIntron, 0.5),1), round(max(SplAnnomat$RatioIntron),1)))

col_fun_prop17 = colorRamp2(c(0, 0.5, 1), c(colwheel[32], "white", colwheel[16]))
lgd17 = Legend(col_fun = col_fun_prop17, title = "PET RIN", at = c(0, 0.5, 1), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
               labels = c(round(min(SplAnnomat$rinvalue),1), round(quantile(SplAnnomat$rinvalue, 0.5),1), round(max(SplAnnomat$rinvalue),1)))

column_ha4 = HeatmapAnnotation(Tissue = SplAnnomat$Tissue,
                               "PET PCT_CODING_BASES" = SplAnnomat$PCT_CODING_BASES_ecdf,
                               "PET PCT_INTRONIC_BASES" = SplAnnomat$PCT_INTRONIC_BASES_ecdf,
                               "PET PCT_INTERGENIC_BASES" = SplAnnomat$PCT_INTERGENIC_BASES_ecdf,
                               "PET PCT_MRNA_BASES" = SplAnnomat$PCT_MRNA_BASES_ecdf,
                               "PET MEDIAN_CV_COVERAGE" = SplAnnomat$MEDIAN_CV_COVERAGE_ecdf,
                               "PET MEDIAN_3PRIME_BIAS" = SplAnnomat$MEDIAN_3PRIME_BIAS_ecdf,
                               "PET MEDIAN_5PRIME_BIAS" = SplAnnomat$MEDIAN_5PRIME_BIAS_ecdf,
                               "PET RatioIntron" = SplAnnomat$RatioIntron_ecdf,
                               "PET RIN" = SplAnnomat$rinvalue_ecdf,
                               annotation_name_gp= gpar(fontsize = 8), annotation_name_rot = c(00),
                               col = list(Tissue = col_ts,
                                          "PET PCT_CODING_BASES" = col_fun_prop9,
                                          "PET PCT_INTRONIC_BASES" = col_fun_prop10,
                                          "PET PCT_INTERGENIC_BASES" = col_fun_prop11,
                                          "PET PCT_MRNA_BASES" = col_fun_prop12,
                                          "PET MEDIAN_CV_COVERAGE" = col_fun_prop13,
                                          "PET MEDIAN_3PRIME_BIAS" = col_fun_prop14,
                                          "PET MEDIAN_5PRIME_BIAS" = col_fun_prop15,
                                          "PET RatioIntron" = col_fun_prop16,
                                          "PET RIN" = col_fun_prop17), show_legend = FALSE,
                               annotation_name_side = c("right"))

# Find two max values (red)
grep(pattern="TCGA-A6-2674",x=colnames(Ratiomat)) # 10
grep(pattern="TCGA-A6-3809",x=colnames(Ratiomat)) # 12

# Find abnormal PET samples (blue)
# TCGA.44.6147 TCGA.A6.5656 TCGA.A6.5659 TCGA.A6.6780 TCGA.A6.6781 TCGA.A7.A13E
# 1.403660     1.403600     1.440131     1.391037     1.418438     1.403342
grep(pattern="TCGA-44-6147",x=colnames(Ratiomat)) # 8
grep(pattern="TCGA-A6-5656",x=colnames(Ratiomat)) # 14
grep(pattern="TCGA-A6-5659",x=colnames(Ratiomat)) # 15
grep(pattern="TCGA-A6-6780",x=colnames(Ratiomat)) # 17
grep(pattern="TCGA-A6-6781",x=colnames(Ratiomat)) # 18
grep(pattern="TCGA-A7-A13E",x=colnames(Ratiomat)) # 21


## ---- 2.3 PET - FFT heatmap

hm40 = Heatmap(Ratiomat, name = "PET - FFT", left_annotation = row_ha5, bottom_annotation = column_ha4,
               clustering_method_rows = "ward.D", show_row_dend = FALSE, row_split = 3,
               clustering_method_columns = "ward.D", column_split = 3,
               show_row_names = FALSE, show_column_names = TRUE,
               use_raster = FALSE,
               column_names_gp = gpar(fontsize=8, col = c(rep("black", 7), c(rep("blue", 1), c(rep("black", 1), rep("red", 1), rep("black", 1), rep("red", 1),
                                                          rep("black", 1), rep("blue", 2), rep("black", 1), rep("blue", 2), rep("black", 2), rep("blue", 1), rep("black", 13))))),
               row_title = "24,798 genes", row_title_gp = gpar(fontface = "bold"),
               column_title = "34 patients \nClustering method: ward.D", column_title_side = "bottom", column_title_gp = gpar(fontface = "bold"),
               show_heatmap_legend = FALSE,
               width = unit(15, "cm"), height = unit(15, "cm"))

# Legend for Ratiomat
col_fun_prop0 = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
lgd0 = Legend(col_fun = col_fun_prop0, title = "PET - FFT", at = c(-4, 0, 4), direction = "horizontal", legend_height = unit(0.5, "mm"), legend_width = unit(25, "mm"), title_position = "lefttop", labels_gp = gpar(fontsize = 8),
              labels = c(round(min(Ratiomat),1), 0, round(max(Ratiomat),1)))


png(file="09_RatioHeatmap_output/RatioHeatmapAnno_ward.D21.png",
    width=1500, height=1000, res = 100)
draw(hm40)

draw(lgd0, x = unit(26.84+1+1+0.01+5, "cm"), y = unit(23.5+1, "cm"))

draw(lgd1, x = unit(26.7+1+1+0.01+5-3, "cm"), y = unit(23-1+1-0.05, "cm")) # Subcategory
draw(lgd18, x = unit(26.7+1+1+0.01+5+1.1, "cm"), y = unit(23-1+1-0.04-0.21, "cm")) # Tissue
draw(lgd2, x = unit(27.89+1+5, "cm"), y = unit(21.4+0.2-1-0.15+1-0.5, "cm"))
draw(lgd3, x = unit(27.87+1+5, "cm"), y = unit(20.9+0.2-1-0.15-0.3+1-0.5, "cm"))
draw(lgd4, x = unit(27.71+1+5, "cm"), y = unit(20.4+0.2-1-0.15-0.6+1-0.5, "cm"))
draw(lgd5, x = unit(27.18+1+5, "cm"), y = unit(19.9+0.2-1-0.15-0.9+1-0.5, "cm"))
draw(lgd6, x = unit(27.18+1+5, "cm"), y = unit(19.4+0.2-1-0.15-1.2+1-0.5, "cm"))
draw(lgd7, x = unit(27+1+5, "cm"), y = unit(18.9+0.2-1-0.15-1.5+1-0.5, "cm"))
draw(lgd8, x = unit(26.93+1+5, "cm"), y = unit(18.4+0.2-1-0.15-1.8+1-0.5, "cm"))

draw(lgd9, x = unit(26.48+1+5, "cm"), y = unit(18.4-0.5+0.2-1-0.15-2.1+1-0.5, "cm"))
draw(lgd10, x = unit(26.34+1+5, "cm"), y = unit(18.4-1+0.2-1-0.15-2.4+1-0.5, "cm"))
draw(lgd11, x = unit(26.1+1+5, "cm"), y = unit(18.4-1.5+0.2-1-0.15-2.7+1-0.5, "cm"))
draw(lgd12, x = unit(26.65+1+5, "cm"), y = unit(18.4-2+0.2-1-0.15-3+1-0.5, "cm"))
draw(lgd13, x = unit(26.21+1+5, "cm"), y = unit(18.4-2.5+0.2-1-0.15-3.3+1-0.5, "cm"))
draw(lgd14, x = unit(26.36+1+5, "cm"), y = unit(18.4-3+0.2-1-0.15-3.6+1-0.5, "cm"))
draw(lgd15, x = unit(26.36+1+5, "cm"), y = unit(18.4-3.5+0.2-1-0.15-3.9+1-0.5, "cm"))
draw(lgd16, x = unit(27.39+1+5, "cm"), y = unit(18.4-4+0.2-1-0.15-4.2+1-0.5, "cm"))
draw(lgd17, x = unit(28.02+1+5, "cm"), y = unit(18.4-4.5+0.2-1-0.15-4.5+1-0.5, "cm"))

dev.off()
