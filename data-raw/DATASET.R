## -------------------------------------------------
##    code to prepare `DATASET` dataset goes here
## -------------------------------------------------

## ---- output

# GeneInfo2 is a 56052 x 22 matrix
# GeneInfo_selected2 is a 101 x 22 matrix
# GeneInfo_KeepTrue2 is a 24798 x 23 matrix

# SampleInfo is a 102 x 39 matrix

# FCmat is a 56052 x 95 matrix
# FCmat_FPKM is a 56052 x 94 matrix
# FCmat_FPKMUQ is a 56052 x 94 matrix
# FCmat_TPM is a 56052 x 94 matrix
# FCmat_TMM_SepByPair_log2 is a 24798 x 94 matrix


## ---- renv

root <- "/Users/myeon/Documents/GitHub/RNAdegrProjR/"
devtools::load_all(root, quiet = TRUE)
devtools::load_all()

## ---- prepare the workspace

# imports
suppressPackageStartupMessages({
  library(dplyr) # For manipulating the data
  library(data.table) # For working with data.tables
  library(biomaRt) # For getting gene information
  # library(org.Hs.eg.db) # For extracting entrezIDs
  library(edgeR) # For TMM normalization
  library(limma) # For TMM normalization
})


## ---- 1. Get gene information table

## ---- 1.1 Gene information table - all 56052 genes

## Gene level data
# Gene length, gene_category
FFvsFFPE <- read.csv(file="data-raw/FFvsFFPE_gene_counts_2nd_strand_total.csv")
dim(FFvsFFPE) # 56052   126

# Add gene_category variable
row_data_add <- read.delim("data-raw/row_data_add.txt")
dim(row_data_add) # 56052     7


## Exon level data
exons <- read.delim("data-raw/gencode.v36.merged.exons.nucleotide.contents.txt")
dim(exons) # 349932     15

# Sum by gene_id
sum_A <- aggregate(exons$X9_num_A, by=list(exons$X4_usercol), FUN=sum)
sum_C <- aggregate(exons$X10_num_C, by=list(exons$X4_usercol), FUN=sum)
sum_G <- aggregate(exons$X11_num_G, by=list(exons$X4_usercol), FUN=sum)
sum_T <- aggregate(exons$X12_num_T, by=list(exons$X4_usercol), FUN=sum)
sum_len <- aggregate(exons$X15_seq_len, by=list(exons$X4_usercol), FUN=sum)

exons1 <- merge(sum_A, sum_C, by="Group.1", all.x=TRUE, sort=TRUE)
exons2 <- merge(exons1, sum_G, by="Group.1", all.x=TRUE, sort=TRUE)
exons3 <- merge(exons2, sum_T, by="Group.1", all.x=TRUE, sort=TRUE)
exons4 <- merge(exons3, sum_len, by="Group.1", all.x=TRUE, sort=TRUE)
colnames(exons4) <- c("gene_id", "sum_A", "sum_C", "sum_G", "sum_T", "sum_len")

# Convert to numeric column in a data frame
exons4 = convert_Chr2Numcol(exons4, 2:6)

# Percent
exons4$exon.pct_A = exons4$sum_A/exons4$sum_len
exons4$exon.pct_C = exons4$sum_C/exons4$sum_len
exons4$exon.pct_G = exons4$sum_G/exons4$sum_len
exons4$exon.pct_T = exons4$sum_T/exons4$sum_len
exons5 <- exons4[,c("gene_id", "exon.pct_A", "exon.pct_C", "exon.pct_G", "exon.pct_T", "sum_len")]

# Weighted sum using seq_len
exons$mX7_pct_at = exons$X7_pct_at*exons$X15_seq_len
exons$mX8_pct_gc = exons$X8_pct_gc*exons$X15_seq_len

wtsum_at <- aggregate(exons$mX7_pct_at, by=list(exons$X4_usercol), FUN=sum)
wtsum_gc <- aggregate(exons$mX8_pct_gc, by=list(exons$X4_usercol), FUN=sum)

exons11 <- merge(wtsum_at, wtsum_gc, by="Group.1", all.x=TRUE, sort=TRUE)
colnames(exons11) <- c("gene_id", "exon.wtsum_at", "exon.wtsum_gc")

# Merge by gene_id
GeneInfo1 <- merge(row_data_add, exons5, by="gene_id", all.x=TRUE, sort=TRUE)
GeneInfo2 <- merge(GeneInfo1, exons11, by="gene_id", all.x=TRUE, sort=TRUE)
GeneInfo2$exon.wtpct_at = GeneInfo2$exon.wtsum_at/GeneInfo2$sum_len
GeneInfo2$exon.wtpct_gc = GeneInfo2$exon.wtsum_gc/GeneInfo2$sum_len
GeneInfo3 <- GeneInfo2[,-c(12:14)] # remove "sum_len", "exon.wtsum_at", "exon.wtsum_gc"
dim(GeneInfo3) # 56052    13


## Intron level data
introns <- read.delim("data-raw/gencode.v36.merged.introns.nucleotide.contents.txt")
dim(introns) # 190140     15

# Sum by gene_id
sum_A <- aggregate(introns$X9_num_A, by=list(introns$X5_usercol), FUN=sum)
sum_C <- aggregate(introns$X10_num_C, by=list(introns$X5_usercol), FUN=sum)
sum_G <- aggregate(introns$X11_num_G, by=list(introns$X5_usercol), FUN=sum)
sum_T <- aggregate(introns$X12_num_T, by=list(introns$X5_usercol), FUN=sum)
sum_len <- aggregate(introns$X15_seq_len, by=list(introns$X5_usercol), FUN=sum)

introns1 <- merge(sum_A, sum_C, by="Group.1", all.x=TRUE, sort=TRUE)
introns2 <- merge(introns1, sum_G, by="Group.1", all.x=TRUE, sort=TRUE)
introns3 <- merge(introns2, sum_T, by="Group.1", all.x=TRUE, sort=TRUE)
introns4 <- merge(introns3, sum_len, by="Group.1", all.x=TRUE, sort=TRUE)
colnames(introns4) <- c("gene_id", "sum_A", "sum_C", "sum_G", "sum_T", "sum_len")

# Convert to numeric column in a data frame
introns4 = convert_Chr2Numcol(introns4, 2:6)

# Percent
introns4$intron.pct_A = introns4$sum_A/introns4$sum_len
introns4$intron.pct_C = introns4$sum_C/introns4$sum_len
introns4$intron.pct_G = introns4$sum_G/introns4$sum_len
introns4$intron.pct_T = introns4$sum_T/introns4$sum_len
introns5 <- introns4[,c("gene_id", "intron.pct_A", "intron.pct_C", "intron.pct_G", "intron.pct_T", "sum_len")]

# Weighted sum using seq_len
introns$mX7_pct_at = introns$X7_pct_at*introns$X15_seq_len
introns$mX8_pct_gc = introns$X8_pct_gc*introns$X15_seq_len

wtsum_at <- aggregate(introns$mX7_pct_at, by=list(introns$X5_usercol), FUN=sum)
wtsum_gc <- aggregate(introns$mX8_pct_gc, by=list(introns$X5_usercol), FUN=sum)

introns11 <- merge(wtsum_at, wtsum_gc, by="Group.1", all.x=TRUE, sort=TRUE)
colnames(introns11) <- c("gene_id", "intron.wtsum_at", "intron.wtsum_gc")

# Merge by gene_id
GeneInfo4 <- merge(GeneInfo3, introns5, by="gene_id", all.x=TRUE, sort=TRUE)
GeneInfo5 <- merge(GeneInfo4, introns11, by="gene_id", all.x=TRUE, sort=TRUE)
GeneInfo5$intron.wtpct_at = GeneInfo5$intron.wtsum_at/GeneInfo5$sum_len
GeneInfo5$intron.wtpct_gc = GeneInfo5$intron.wtsum_gc/GeneInfo5$sum_len
GeneInfo6 <- GeneInfo5[,-c(18:20)] # remove "sum_len", "intron.wtsum_at", "intron.wtsum_gc"
dim(GeneInfo6) # 56052    19


## Transcript level data
utr <- read.delim("data-raw/gencode.v36.utr.txt", header=FALSE)
dim(utr) # 169127      7

# Length by gene_id_UTR
# length = |end point-start point| + 1
utr$len = abs(utr$V3-utr$V2)+1
utr$gene_id_UTR = paste(utr$V5, utr$V6, sep = "_")

utr1 <- aggregate(utr$len, by=list(utr$gene_id_UTR), FUN=mean)
utr1$gene_id <- sub("_.*", "", utr1$Group.1)
utr1$UTR <- sub(".*_", "", utr1$Group.1)
colnames(utr1) <- c("Group.1", "meanlen", "gene_id", "UTR")

utr21 <- utr1[utr1$UTR=="3UTR",c("gene_id","meanlen")]
colnames(utr21) <- c("gene_id","3UTR.meanlen")
utr22 <- utr1[utr1$UTR=="5UTR",c("gene_id","meanlen")]
colnames(utr22) <- c("gene_id","5UTR.meanlen")
utr3 <- merge(utr21, utr22, by="gene_id", all.x=TRUE, sort=TRUE)
dim(utr3) # 19977     3

# Remove duplicates
utr4 <- utr3[!duplicated(utr3$gene_id), ]
dim(utr4) # 19925     3

# Add subgroup
subgroup <- read.delim("data-raw/GeneInfo_subgroup_annotation.txt")

# Merge by gene_id
GeneInfo1 <- merge(GeneInfo6, utr4, by="gene_id", all.x=TRUE, sort=TRUE)
GeneInfo2 <- merge(GeneInfo1, subgroup, by="gene_category", all.x=TRUE, sort=TRUE)
dim(GeneInfo2) # 56052    22


## ---- 1.2 Gene information table - a subset of 101 genes

# Keep 101 genes
Allgenes <- list.files("data-raw/pileup_selected")
genes <- sub("_rs_pileup_part_intron.RData", "", Allgenes)
length(genes) # 101

# Find exact match
keeprownum = matrix(0, length(genes), 1)
for (k in 1:length(genes)){
  kthgene <- genes[k]
  keeprownum[k,] = grep(pattern=sprintf("^%s$",kthgene),x=c(GeneInfo2[,3]))
}
GeneInfo_selected2 <- GeneInfo2[keeprownum,]
dim(GeneInfo_selected2) # 101  22


## ---- 1.3 Gene information table - a subset of 24798 genes

load("data-raw/Keep.RData")
dim(keep1) # 56052     1

keep2 <- cbind(rownames(keep1), keep1)
GeneInfo_KeepTrue1 <- merge(GeneInfo2, keep2, by.x="geneSymbol", by.y="V1", all.x=TRUE, sort=TRUE)

# Make the selected variable
keeprownum = matrix(0, length(genes), 1)
for (k in 1:length(genes)){
  kthgene <- genes[k]
  keeprownum[k,] = grep(pattern=sprintf("^%s$",kthgene),x=c(GeneInfo_KeepTrue1[,1]))
}
GeneInfo_KeepTrue1$selected <- 0
GeneInfo_KeepTrue1[keeprownum,"selected"] <- 1

# Keep 24798 genes
GeneInfo_KeepTrue2 <- GeneInfo_KeepTrue1[GeneInfo_KeepTrue1$V2=="TRUE",]
GeneInfo_KeepTrue2 <- GeneInfo_KeepTrue2[,-c(23)]
dim(GeneInfo_KeepTrue2) # 24798    23


## ---- 2. Get sample information table

## Pilot2_entityID
Pilot <- read.csv(file="data-raw/base/Pilot2_entityID.csv")
dim(Pilot) # 120   8

Pilot$part1 <- substr(Pilot$entity_submitter_id, start=1, stop=12)
Pilot$part2[Pilot$Preservation=="Fresh Frozen" & Pilot$Protocol=="mRNA"] <- "FFM"
Pilot$part2[Pilot$Preservation=="Fresh Frozen" & Pilot$Protocol=="TotalRNAseq"] <- "FFT"
Pilot$part2[Pilot$Preservation=="FFPE" & Pilot$Protocol=="TotalRNAseq"] <- "PET"
Pilot$NewSampleId <- paste(Pilot$part1, Pilot$part2, sep = "-")

## Select samples with all 3 pairs -> make a new file with these samples
# Keep Primary TissueType
Pilot0 <- Pilot[Pilot$TissueType=="Primary", ]

# Keep all 3 pairs for each patient
numsample <- aggregate(Pilot0$NewSampleId, by=list(Pilot0$part1), FUN=length)
Pilot1 <- merge(Pilot0, numsample, by.x="part1", by.y="Group.1")
Pilot2 <- Pilot1[Pilot1$x=="3", ]
dim(Pilot2) # 102  12

# Drop unnecessary columns
Pilot3 <- Pilot2[,-c(1,10,12)]
dim(Pilot3) # 102   9


## Picard tool CollectRnaSeqMetrics, METRICS CLASS
Allfiles <- list.files("data-raw/FFvsFFPE_ream_matrix")
files <- sub(".RNA_Metrics", "", Allfiles)

f=1
fthfile <- files[f]
all_content = readLines(sprintf("data-raw/FFvsFFPE_ream_matrix/%s.RNA_Metrics", fthfile))
keep_rows = all_content[c(7:8)]
matclass <- read.table(textConnection(keep_rows), sep="\t", header=T, na.strings="", stringsAsFactors=F)
matclass1 <- t(matclass)

Mmat = matrix(0, nrow(matclass1),length(files)) # a 30 x 120 matrix
for (f in 1:length(files)){
  fthfile <- files[f]
  all_content = readLines(sprintf("data-raw/FFvsFFPE_ream_matrix/%s.RNA_Metrics", fthfile))
  keep_rows = all_content[c(7:8)]
  matclass <- read.table(textConnection(keep_rows), sep="\t", header=T, na.strings="", stringsAsFactors=F)
  matclass1 <- t(matclass)

  Mmat[,f] <- as.matrix(matclass1)
}
colnames(Mmat) <- c(files)
rownames(Mmat) <- rownames(matclass1)

# Keep the selected patients in Pilot data
MCmat = rpl_NewSampleId(mat=Mmat, Pilot2=Pilot2)
dim(MCmat) # 30 102


## Degradation indexes: RatioIntron = INTRONIC_BASES/CODING_BASES
MCmat1 <- t(MCmat[-c(28:30),])
MCmat2 <- as.matrix(MCmat1[,6]/MCmat1[,4])

# Merge RatioIntron
MCmat11 <- cbind(rownames(MCmat1), MCmat1)
MCmat21 <- cbind(rownames(MCmat2), MCmat2)
MCmat3 <- merge(MCmat11, MCmat21, by="V1", all.x=TRUE, sort=TRUE)
colnames(MCmat3) <- c("NewSampleId", colnames(MCmat3[,-c(1,ncol(MCmat3))]), "RatioIntron")
dim(MCmat3) # 102  29


## Degradation indexes: ratio_28s_18s, rinvalue
RatioRIN <- read.csv(file="data-raw/Pilot2_entityID_filtered_with_RIN.csv")
dim(RatioRIN) # 102  10

RatioRIN$part1 <- substr(RatioRIN$entity_submitter_id, start=1, stop=12)
RatioRIN$part2[RatioRIN$Preservation=="Fresh Frozen" & RatioRIN$Protocol=="mRNA"] <- "FFM"
RatioRIN$part2[RatioRIN$Preservation=="Fresh Frozen" & RatioRIN$Protocol=="TotalRNAseq"] <- "FFT"
RatioRIN$part2[RatioRIN$Preservation=="FFPE" & RatioRIN$Protocol=="TotalRNAseq"] <- "PET"
RatioRIN$NewSampleId <- paste(RatioRIN$part1, RatioRIN$part2, sep = "-")

RatioRIN1 <- RatioRIN[,c("NewSampleId", "ratio_28s_18s", "rinvalue")]
dim(RatioRIN1) # 102   3

# Merge by NewSampleId
SI1 <- merge(Pilot3, MCmat3, by="NewSampleId", all.x=TRUE, sort=TRUE)
SI2 <- merge(SI1, RatioRIN1, by="NewSampleId", all.x=TRUE, sort=TRUE)

SampleInfo = convert_Chr2Numcol(SI2, 10:37)
dim(SampleInfo) # 102  39


## ---- 3. Gene expression normalization

## ---- 3.1 Gene fragment counts matrix - raw

Allfiles <- list.files("data-raw/countstsv")
files <- sub(".rna_seq.augmented_star_gene_counts.tsv", "", Allfiles)

f=1
fthfile <- files[f]
all_content = readLines(sprintf("data-raw/countstsv/%s.rna_seq.augmented_star_gene_counts.tsv", fthfile))
skip_rows = all_content[-c(1,3:6)]
gene_counts <- read.table(textConnection(skip_rows), sep="\t", header=T, na.strings="", stringsAsFactors=F)

Fmat = matrix(0, dim(gene_counts)[1],length(files)) # dim(gene_counts)[1]=60660
col2name = matrix(0, 1,length(files))
for (f in 1:length(files)){
  fthfile <- files[f]
  all_content = readLines(sprintf("data-raw/countstsv/%s.rna_seq.augmented_star_gene_counts.tsv", fthfile))
  skip_rows = all_content[-c(1,3:6)]
  gene_counts <- read.table(textConnection(skip_rows), sep="\t", header=T, na.strings="", stringsAsFactors=F)

  Fmat[,f] <- as.matrix(gene_counts[ ,"stranded_second"])

  # Find NewSampleId using .tsv file name (same with file_name in Pilot)
  col2name[,f] <- Pilot[which(grepl(pattern=sprintf("%s", fthfile), Pilot$file_name)), "entity_submitter_id"]
}
rownames(Fmat) <- gene_counts[,"gene_id"]
colnames(Fmat) <- c(col2name)

# Keep the selected patients in Pilot data
Fmat0 = rpl_NewSampleId(mat=Fmat, Pilot2=Pilot2) # Rows: gene_id

# Remove 8 NewSampleIds that had no expression in the pileup FFM (FF mRNA-seq)
removeids <- c("TCGA-A6-2674-FFM",
               "TCGA-A6-2684-FFM",
               "TCGA-A6-3809-FFM",
               "TCGA-A6-3810-FFM",
               "TCGA-BK-A0CA-FFM",
               "TCGA-BK-A0CC-FFM",
               "TCGA-BK-A139-FFM",
               "TCGA-BK-A26L-FFM")
removecolnum = matrix(0, 1, 8)
for (r in 1:8){
  removecolnum[,r] = which(grepl(pattern=removeids[r],x=colnames(Fmat0)))
}
Fmat1 <- Fmat0[,-removecolnum]
dim(Fmat0) # 60660   102
dim(Fmat1) # 60660    94(102-8)


## Change row names of counts matrix to have geneSymbol & keep only 56052 genes
# Check if geneSymbol has unique values
FFvsFFPE <- read.csv(file="data-raw/counts/FFvsFFPE_gene_counts_2nd_strand_total.csv")
dim(FFvsFFPE) # 56052   126

AllgeneSymbolLists <- FFvsFFPE$geneSymbol
UnqgeneSymbolLists <- unique(FFvsFFPE$geneSymbol)
length(AllgeneSymbolLists) # 56052
length(UnqgeneSymbolLists) # 56052

# Check if gene_id has unique values
Allgene_idLists <- FFvsFFPE$gene_id
Unqgene_idLists <- unique(FFvsFFPE$gene_id)
length(Allgene_idLists) # 56052
length(Unqgene_idLists) # 56052

# Keep only 56052 genes
Fmat2 <- cbind(rownames(Fmat1), Fmat1)
Fmat3 <- merge(FFvsFFPE[,c("gene_id","geneSymbol","merged")], Fmat2, by.x="gene_id", by.y="V1", all.x=TRUE, sort=TRUE)

# Change the unit of gene length to kb
merged_kb = Fmat3$merged/1000
Fmat4 <- cbind(merged_kb, Fmat3)

# Sort by geneSymbol
Fmat5 <- Fmat4[order(Fmat4$geneSymbol),]

# Change rows from gene_id to geneSymbol
Fmat6 <- Fmat5[,-c(2:4)] # remove "gene_id","geneSymbol","merged"

# Convert dataframe column from character to numeric
FCmat = convert_Chr2Numcol(Fmat6, 2:ncol(Fmat6))
rownames(FCmat) <- c(Fmat5$geneSymbol)
dim(FCmat) # 56052    95


## ---- 3.2 Gene fragment counts matrix - FPKM normalization

ScalingFactor <- matrix(apply(FCmat[,-1],2,sum)/10^6, 1, ncol(FCmat[,-1]))

SFmat <- repmat(ScalingFactor, nrow(FCmat), 1)
dim(SFmat) # 56052    94
FCmat1 <- FCmat[,-1]/SFmat
FCmat_FPKM <- FCmat1/repmat(matrix(merged_kb,nrow(FCmat),1),1,ncol(FCmat1))
dim(FCmat_FPKM) # 56052    94


## ---- 3.3 Gene fragment counts matrix - FPKM-UQ normalization

# Change sum(C_i) to C_qtl(0.75)*G where G is nrow(FCmat)=56052
ScalingFactor <- matrix((apply(FCmat[,-1],2,quantile)[4,]*nrow(FCmat))/10^6, 1, ncol(FCmat[,-1]))

SFmat <- repmat(ScalingFactor, nrow(FCmat), 1)
dim(SFmat) # 56052    94
FCmat2 <- FCmat[,-1]/SFmat
FCmat_FPKMUQ <- FCmat2/repmat(matrix(merged_kb,nrow(FCmat),1),1,ncol(FCmat2))
dim(FCmat_FPKMUQ) # 56052    94


## ---- 3.4 Gene fragment counts matrix - TPM normalization

deno <- repmat(matrix(apply(FCmat_FPKM,2,sum),1,ncol(FCmat_FPKM)), nrow(FCmat_FPKM), 1)
FCmat_TPM <- (FCmat_FPKM/deno)*10^6
dim(FCmat_TPM) # 56052    94


## ---- 3.5 Gene fragment counts matrix - TMM normalization

# Design matrix
merged <- FCmat[,"merged_kb"]*1000
FCmat1 <- cbind(merged, FCmat[,-1])
pair <- substr(colnames(FCmat1[,-1]), start=14, stop=16)

pair <- factor(pair)
table(pair)
# pair
# FFM FFT PET
# 26  34  34

design <- model.matrix(~0+pair)
colnames(design) <- levels(pair)
design

# Filtering to remove low counts
y <- FCmat1
keep <- filterByExpr(y, design, min.count=10, min.prop=0.7)
table(keep)
# keep
# FALSE  TRUE
# 31254 24798

y <- y[keep,]
dim(y) # 24798    95

# Seperate a matrix by pair, then nomalize TMM
pair2 <- substr(colnames(y), start=14, stop=16)
keepids <- c("FFM","FFT","PET")

FFMids <- which(grepl(pattern=keepids[1],x=pair2))
FFTids <- which(grepl(pattern=keepids[2],x=pair2))
PETids <- which(grepl(pattern=keepids[3],x=pair2))

FFMmat <- y[,c(1,FFMids)]
FFTmat <- y[,c(1,FFTids)]
PETmat <- y[,c(1,PETids)]

# Calculate TMM; convert the read counts into log2(x+1)
y1 <- DGEList(FFMmat[,-1], genes=FFMmat[,1,drop=FALSE])
options(digits=3)
y1 <- calcNormFactors(y1)
FCmat_TMM_SepByPair1_log2 <- log2(cpm(y1, method="TMM", refColumn=1)+1)

y2 <- DGEList(FFTmat[,-1], genes=FFTmat[,1,drop=FALSE])
options(digits=3)
y2 <- calcNormFactors(y2)
FCmat_TMM_SepByPair2_log2 <- log2(cpm(y2, method="TMM", refColumn=1)+1)

y3 <- DGEList(PETmat[,-1], genes=PETmat[,1,drop=FALSE])
options(digits=3)
y3 <- calcNormFactors(y3)
FCmat_TMM_SepByPair3_log2 <- log2(cpm(y3, method="TMM", refColumn=1)+1)

FCmat_TMM_SepByPair_log2 <- cbind(FCmat_TMM_SepByPair1_log2, FCmat_TMM_SepByPair2_log2, FCmat_TMM_SepByPair3_log2)
dim(FCmat_TMM_SepByPair_log2) # 24798    94


## --------------------------------------
##     save processed data to data/
## --------------------------------------

usethis::use_data(GeneInfo2, overwrite = TRUE)   # gene information table - all 56052 genes
usethis::use_data(GeneInfo_selected2, overwrite = TRUE)   # gene information table - a subset of 101 genes
usethis::use_data(GeneInfo_KeepTrue2, overwrite = TRUE)   # gene information table - a subset of 24798 genes

usethis::use_data(SampleInfo, overwrite = TRUE)   # sample information table

usethis::use_data(FCmat, overwrite = TRUE)   # gene fragment counts matrix - raw
usethis::use_data(FCmat_FPKM, overwrite = TRUE)   # gene fragment counts matrix - FPKM normalization
usethis::use_data(FCmat_FPKMUQ, overwrite = TRUE)   # gene fragment counts matrix - FPKM-UQ normalization
usethis::use_data(FCmat_TPM, overwrite = TRUE)   # gene fragment counts matrix - TPM normalization
usethis::use_data(FCmat_TMM_SepByPair_log2, overwrite = TRUE)   # gene fragment counts matrix - TMM normalization
