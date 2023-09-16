## ---------------------------------------------
##    code for GC percentage TMM scatter plot
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
  library(ggplot2) # For ggplot
  library(gridExtra) # For multiple ggplots in pdf
})


## ---- 1. TMM matrix

load("data/FCmat_TMM_SepByPair_log2.rda")
dim(FCmat_TMM_SepByPair_log2) # 24798    94

# Separate columns by pair
pair <- substr(colnames(FCmat_TMM_SepByPair_log2), start=14, stop=16)
keepcol <- c("FFM","FFT","PET")

FFTcol <- grep(pattern=keepcol[2],x=pair)
PETcol <- grep(pattern=keepcol[3],x=pair)

FFTmat <- FCmat_TMM_SepByPair_log2[,c(FFTcol)]
PETmat <- FCmat_TMM_SepByPair_log2[,c(PETcol)]
dim(FFTmat) # 24798    34
dim(PETmat) # 24798    34

# Sort by colnames (patients)
FFTmat1 <- FFTmat[, order(colnames(FFTmat))]
PETmat1 <- PETmat[, order(colnames(PETmat))]

# Sort by rownames (genes)
FFTmat2 <- FFTmat1[order(rownames(FFTmat1)), ]
PETmat2 <- PETmat1[order(rownames(PETmat1)), ]
dim(FFTmat2) # 24798    34
dim(PETmat2) # 24798    34


## Filter by mad > quantile (keep top 75% base on mad in FFT)
vec.mad = as.matrix(apply(FFTmat2, 1, mad))
dim(vec.mad) # 24798     1

vec.mad1 = as.matrix(vec.mad[vec.mad>quantile(vec.mad, probs = c(.75)), ])
dim(vec.mad1) # 6200    1

FFTmat3 <- FFTmat2[rownames(vec.mad1),]
PETmat3 <- PETmat2[rownames(vec.mad1),]
dim(FFTmat3) # 6200   34
dim(PETmat3) # 6200   34

FFTmat4 = convert_pivot.longer(FFTmat3, c("gene", "sample", "TMM"))
dim(FFTmat4) # 210800      3
PETmat4 = convert_pivot.longer(PETmat3, c("gene", "sample", "TMM"))
dim(PETmat4) # 210800      3
dim(FFTmat3)[1]*dim(FFTmat3)[2] # 210800


## Add gene info
load("data/GeneInfo_KeepTrue2.rda")
dim(GeneInfo_KeepTrue2) # 24798    23

pair <- c(rep("FFT", dim(FFTmat3)[1]*dim(FFTmat3)[2]), rep("PET", dim(FFTmat3)[1]*dim(FFTmat3)[2]))
mat <- cbind(pair, rbind(FFTmat4, PETmat4))

mat2 <- merge(mat, GeneInfo_KeepTrue2[,c("geneSymbol","exon.wtpct_gc","intron.wtpct_gc")], by.x="gene", by.y="geneSymbol", all.x=TRUE, sort=TRUE)
dim(mat2) # 421600      6

# Convert . to - in sample for NewSampleId
mat2$NewSampleId <- gsub("[.]", "-", mat2$sample)
mat2$id <- substr(mat2$NewSampleId, start=1, stop=12)


## ---- 2. Scatter plot: All subcategory

## Sample info for sorting by RatioIntron in PET
load("data/SampleInfo.rda")
SampleInfo <- SampleInfo[substr(SampleInfo$NewSampleId, start=14, stop=16)=="PET", ]
dim(SampleInfo) # 34 39

SampleInfo1 <- SampleInfo[,c("NewSampleId","RatioIntron")]
dim(SampleInfo1) # 34  2

# Start from the max value of RatioIntron
SampleInfo2 <- SampleInfo1[order(SampleInfo1$RatioIntron, decreasing=TRUE), ]

Allpt <- unique(substr(SampleInfo2$NewSampleId, start=1, stop=12))
length(Allpt) # 34

RI <- round(SampleInfo2$RatioIntron,2)

P1 <- list()
P2 <- list()
for(i in 1:length(Allpt)) {
  P1[[i]] <- list()
  P2[[i]] <- list()

  P1[[i]] <- ggplot(mat2[mat2$id==Allpt[i], ], aes(x=exon.wtpct_gc, y=TMM, color=pair)) +
    geom_point(aes(y=TMM, alpha=0.05)) +
    scale_color_manual(values = c("FFT" = "blue", "PET" = "green")) +
    labs(title=paste0(Allpt[i], ", RatioIntron=",RI[i]), x="Exon GC percentage", y="log2(TMM+1)") +
    theme(legend.position="bottom") +
    theme(panel.background = element_rect(fill="gray97")) +
    geom_smooth(method="loess", formula=y~x, se=F)

  P2[[i]] <- ggplot(mat2[mat2$id==Allpt[i], ], aes(x=intron.wtpct_gc, y=TMM, color=pair)) +
    geom_point(aes(y=TMM, alpha=0.05)) +
    scale_color_manual(values = c("FFT" = "blue", "PET" = "green")) +
    labs(title=paste0(Allpt[i], ", RatioIntron=",RI[i]), x="Intron GC percentage", y="log2(TMM+1)") +
    theme(legend.position="bottom") +
    theme(panel.background = element_rect(fill="gray97")) +
    geom_smooth(method="loess", formula=y~x, se=F)
}

pdf(file="09_RatioHeatmap_output/GCpctTMMScatterplot_AllSubcategoryMad75.pdf", width=12, height=9,
    bg="white", colormodel="cmyk", paper="letter", onefile=TRUE)
grid.arrange(grobs=list(P1[[1]],P2[[1]], P1[[2]],P2[[2]]), nrow=2)
grid.arrange(grobs=list(P1[[3]],P2[[3]], P1[[4]],P2[[4]]), nrow=2)
grid.arrange(grobs=list(P1[[5]],P2[[5]], P1[[6]],P2[[6]]), nrow=2)
grid.arrange(grobs=list(P1[[7]],P2[[7]], P1[[8]],P2[[8]]), nrow=2)
grid.arrange(grobs=list(P1[[9]],P2[[9]], P1[[10]],P2[[10]]), nrow=2)

grid.arrange(grobs=list(P1[[11]],P2[[11]], P1[[12]],P2[[12]]), nrow=2)
grid.arrange(grobs=list(P1[[13]],P2[[13]], P1[[14]],P2[[14]]), nrow=2)
grid.arrange(grobs=list(P1[[15]],P2[[15]], P1[[16]],P2[[16]]), nrow=2)
grid.arrange(grobs=list(P1[[17]],P2[[17]], P1[[18]],P2[[18]]), nrow=2)
grid.arrange(grobs=list(P1[[19]],P2[[19]], P1[[20]],P2[[20]]), nrow=2)

grid.arrange(grobs=list(P1[[21]],P2[[21]], P1[[22]],P2[[22]]), nrow=2)
grid.arrange(grobs=list(P1[[23]],P2[[23]], P1[[24]],P2[[24]]), nrow=2)
grid.arrange(grobs=list(P1[[25]],P2[[25]], P1[[26]],P2[[26]]), nrow=2)
grid.arrange(grobs=list(P1[[27]],P2[[27]], P1[[28]],P2[[28]]), nrow=2)
grid.arrange(grobs=list(P1[[29]],P2[[29]], P1[[30]],P2[[30]]), nrow=2)

grid.arrange(grobs=list(P1[[31]],P2[[31]], P1[[32]],P2[[32]]), nrow=2)
grid.arrange(grobs=list(P1[[33]],P2[[33]], P1[[34]],P2[[34]]), nrow=2)

dev.off()
