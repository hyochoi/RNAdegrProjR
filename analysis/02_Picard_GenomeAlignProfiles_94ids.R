## ----------------------------------------
##    code for genome alignment profiles
## ----------------------------------------

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


## ---- 1. Outputs from Picard tool CollectRnaSeqMetrics, METRICS CLASS

load("data-raw/Pilot2.RData") # from DATASET.R

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


## ---- 2. Genome alignment profiles for 94 samples

# Remove FFM 8 samples that had no expression in the pileup FFM (FF mRNA-seq)
removespls <- c("TCGA-A6-2674-FFM",
                "TCGA-A6-2684-FFM",
                "TCGA-A6-3809-FFM",
                "TCGA-A6-3810-FFM",
                "TCGA-BK-A0CA-FFM",
                "TCGA-BK-A0CC-FFM",
                "TCGA-BK-A139-FFM",
                "TCGA-BK-A26L-FFM")
removecolnum = matrix(0, 1, 8)
for (r in 1:8){
  removecolnum[,r] = which(grepl(pattern=removespls[r],x=colnames(MCmat)))
}
MCmat0 <- MCmat[,-removecolnum]
dim(MCmat0) # 30 94(102-8)

MCmat1 <- t(MCmat0[1:7,])
Unaligned = MCmat1[,1]-MCmat1[,2] # Unaligned = PF_BASES-PF_ALIGNED_BASES
MCmat2 <- cbind(Unaligned, MCmat1)
MCmat3 <- MCmat2[,-c(2:3)]
SUM <- apply(MCmat3, 1, sum) # SUM = sum(Unaligned, RIBOSOMAL_BASES, CODING_BASES, UTR_BASES, INTRONIC_BASES, INTERGENIC_BASES)
MCmat4 <- cbind(MCmat3, SUM)
colnames(MCmat4)
#[1] "Unaligned"        "RIBOSOMAL_BASES"  "CODING_BASES"     "UTR_BASES"        "INTRONIC_BASES"   "INTERGENIC_BASES" "SUM"

PCT_Unaligned = (MCmat4[,1]/MCmat4[,7])*100
PCT_Intergenic = (MCmat4[,6]/MCmat4[,7])*100
PCT_Intronic = (MCmat4[,5]/MCmat4[,7])*100
PCT_Coding.UTR = ((MCmat4[,3]+MCmat4[,4])/MCmat4[,7])*100
PCT_Ribo = (MCmat4[,2]/MCmat4[,7])*100
PCTmat <- cbind(PCT_Unaligned, PCT_Intergenic, PCT_Intronic, PCT_Coding.UTR, PCT_Ribo)
dim(PCTmat) # 94  5

pair <- substr(rownames(MCmat4), start=14, stop=16)
df <- data.frame(pair, PCTmat)
dim(df) # 94  6


pair1 <- rep(df[,c("pair")],5)
grp1 <- c(rep("Unaligned",nrow(df)),rep("Intergenic",nrow(df)),
        rep("Intronic",nrow(df)), rep("Coding+UTR",nrow(df)),
        rep("Ribo",nrow(df)))
PCT1 <- matrix(c(df[,2],df[,3],df[,4],df[,5],df[,6]), nrow(df)*5, 1)
df1 <- data.frame(pair1, grp1, PCT1)
dim(df1) # 470(94*5)   3

# Remove PCT_Ribo since it closes to 0
MCdf <- df1[grp1!="Ribo",]
dim(MCdf) # 376(94*4)   3


## Violin plot
MCdf$grp1 <- factor(MCdf$grp1, levels=c("Unaligned", "Intergenic", "Intronic", "Coding+UTR"), ordered=TRUE)

p3 <- ggplot(MCdf, aes(x=grp1, y=PCT1, color=pair1)) +
  scale_colour_brewer(palette="Set1") +
  geom_violin(position = position_dodge(width = 0.9)) +
  geom_point(position = position_jitterdodge(seed=1, dodge.width=0.9), alpha=0.4, shape=16, size=1) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, position=position_dodge(width = 0.9)) +
  theme(legend.position="bottom") +
  labs(title="Violin plot by grp",x="", y="") +
  scale_x_discrete(limits = c("Unaligned", "Intergenic", "Intronic", "Coding+UTR")) +
  theme(panel.background = element_rect(fill="gray97")) +
  guides(colour=guide_legend(title="",nrow=1,byrow=TRUE))

p4 <- ggplot(MCdf, aes(x=pair1, y=PCT1, color=grp1)) +
  scale_colour_brewer(palette="Dark2") +
  geom_violin(position = position_dodge(width = 0.9)) +
  geom_point(position = position_jitterdodge(seed=1, dodge.width=0.9), alpha=0.4, shape=16, size=1) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, position=position_dodge(width = 0.9)) +
  theme(legend.position="bottom") +
  labs(title="Violin plot by pair",x="", y="") +
  theme(panel.background = element_rect(fill="gray97")) +
  guides(colour=guide_legend(title="",nrow=1,byrow=TRUE))

pdf(file="02_Picard_output/Picard_GenomeAlignProfiles_94ids.pdf", width=12, height=9,
    bg="white", colormodel="cmyk", paper="letter", onefile=TRUE)
grid.arrange(grobs=list(p3,p4), nrow=2)
dev.off()
