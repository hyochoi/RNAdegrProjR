Demo
================

## Data Description

### Inputs
We need datasets such as `genelist`, coverage `pileup`, and `sampleInfo` to obtain the sample quality index outputs and plots.
`geneInfo` is optional if you want to compare results by the properties of genes. The description of Input data and variable names is listed as follows:

- `genelist`: a vector of gene names
- `pileupPath`: a vector for file paths of coverage pileupData including .RData file names
- `geneInfo`: a data frame of gene information including gene ID and properties based on _gencode v36_
  - `gene_id`: Ensembl gene ID
  - `geneSymbol`: gene names
  - `merged`: gene length
  - `exon.wtpct_gc`: weighted percentage of GC from exon level data
  - `subcategory`: protein-coding or lncRNA
- `sampleInfo`: a data frame of sample information including sample ID and properties from [_Picard RnaSeqMetrics_](https://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics)
  - `SampleID`: sample ID
  - `PF_BASES`: the total number of bases within the PF_READS of the SAM or BAM file to be examined
  - `PF_ALIGNED_BASES`: the total number of aligned bases, in all mapped PF reads, that are aligned to the reference sequence
  - `RIBOSOMAL_BASES`: number of bases in primary alignments that align to ribosomal sequence
  - `CODING_BASES`: number of bases in primary alignments that align to a non-UTR coding base for some gene, and not ribosomal sequence
  - `UTR_BASES`: number of bases in primary alignments that align to a UTR base for some gene, and not a coding base
  - `INTRONIC_BASES`: number of bases in primary alignments that align to an intronic base for some gene, and not a coding or UTR base
  - `INTERGENIC_BASES`: number of bases in primary alignments that do not align to any gene
  - `RINs`: RIN value
  - `RatioIntron`: ratio of intronic bases and coding bases

### Read BAM

### Alliance
This example consists of 1,000 selected genes among protein-coding and lncRNA genes and FFT 171 samples, which can be found in [data](https://github.com/hyochoi/RNAdegrProjR/tree/main/data).

## Genome Alignment Profiles

The relative coverage of reads mapped in unaligned (unmapped bases),
intergenic, intronic, and exonnic/protein coding and UTR regions

``` r
GAP = plot_GAP(sampleInfo, plot=TRUE)
round(apply(GAP$PCTmat, 2, mean), 1)
 # PCT_Unaligned PCT_Intergenic   PCT_Intronic PCT_Coding.UTR 
 #          10.8           20.8           32.9           35.5

print(GAP$plot)
```
![](figures/Allianceex_GAP_v058.png)<!-- -->
<div align="center">
  <img width="40%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/figures/Allianceex_GAP_v058.png">
</div>

## Gene Body Coverage

### Gene body coverage with all samples
The union transcript is used to extract only exon pileup. To keep only exon location, we first build coverage `pileup` from raw pileup (part_intron) to `pileupData` (only_exon). 
Let’s compare the dimension of `pileup` of the first and the last genes using `get_pileupExon` function: _LINC01772_ and _MIR133A1HG_ have 3,245 and 5,825 positions, respectively.

``` r
genelist[c(1, length(genelist))]
```

    ## [1] "LINC01772"  "MIR133A1HG"

``` r
LI = get_pileupExon(g=1, pileupPath)
AC = get_pileupExon(g=length(genelist), pileupPath)
dim(LI); dim(AC)
```

    ## [1] 3245  171

    ## [1] 5825  171

An equal number of positions will be selected per genes, and the evenly spaced regions and pileup are defined as gene body percentile and normalized coverage in the plot in `plot_GBC`.
See R function in Documentation for details. Mean scaled normalized
coverage at the gene body percentile by gene length.

``` r
## Filtered genes
TPM <- read.table(paste0(root2, "Data/Alliance_TPM_proteincoding.txt"),sep="\t",header=T,check.names=F,stringsAsFactors=F)
genelist2 <- filter_lowExpGenes(genelist, TPM, thru=5, pct=50)
pileupPath2 <- paste0(folder_path, "/", genelist2, "/", genelist2,"_pileup_part_intron.RData")
length(pileupPath2) # 788

geneInfo2 <- geneInfo[match(genelist2, geneInfo$geneSymbol), ] %>%
  mutate(merged_kb=merged/1000,
         Len=factor(case_when(
           merged_kb>=0 & merged_kb<5 ~ "0~5 kb",
           merged_kb>=5 ~ "5+ kb",
           is.na(merged_kb) ~ NA)),
         LenSorted=forcats::fct_relevel(Len, "0~5 kb", "5+ kb")
  )

table(geneInfo2$LenSorted)
# 0~5 kb  5+ kb
# 64    724

sum(table(geneInfo2$LenSorted)) # 788

genelist2Len0 <- geneInfo2[geneInfo2$LenSorted=="0~5 kb", c("geneSymbol")]
pileupPath2Len0 <- paste0(folder_path, "/", genelist2Len0, "/", genelist2Len0,"_pileup_part_intron.RData")
length(pileupPath2Len0) # 64

genelist2Len5 <- geneInfo2[geneInfo2$LenSorted=="5+ kb", c("geneSymbol")]
pileupPath2Len5 <- paste0(folder_path, "/", genelist2Len5, "/", genelist2Len5,"_pileup_part_intron.RData")
length(pileupPath2Len5) # 724


GBC0 = plot_GBC(pileupPath2Len0, geneNames=genelist2Len0, rnum=100, method=1, scale=TRUE, stat=2, plot=TRUE, sampleInfo)
GBC5 = plot_GBC(pileupPath2Len5, geneNames=genelist2Len5, rnum=100, method=1, scale=TRUE, stat=2, plot=TRUE, sampleInfo)

p0 <- GBC0$plot +
  coord_cartesian(ylim=c(0, 0.017)) +
  ggtitle("0~5 kb")
p5 <- GBC5$plot +
  coord_cartesian(ylim=c(0, 0.017)) +
  ggtitle("5+ kb")

ggpubr::ggarrange(p0, p5, common.legend=TRUE, legend="bottom", nrow=1)
```
![](figures/Allianceex_GBC_v058.png)<!-- -->
<div align="center">
  <img width="70%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/figures/Allianceex_GBC_v058.png">
</div>

### CVs from gene body coverage

Metrics from scaled normalized transcript coverage for samples

``` r
ptm=proc.time()[3]
met = get_metrics(pileupPath, geneNames=genelist, rnum=100, method=1, scale=TRUE, margin=1)
proc.time()[3]-ptm
#  elapsed 
# 1939.916 

head(met)
#                         mean          sd        CV     median         mad
# S000004-37958-003 0.01034679 0.003265439 0.3155994 0.01054433 0.003117695
# S000004-38071-002 0.01035084 0.003249576 0.3139432 0.01055158 0.003045546
# S000004-38094-003 0.01031615 0.003287909 0.3187146 0.01051358 0.003152118
# S000004-38120-002 0.01048696 0.003516319 0.3353041 0.01065609 0.003372997
# S000004-38134-002 0.01046729 0.003521123 0.3363929 0.01056110 0.003305320
# S000004-38138-003 0.01067134 0.003933500 0.3686042 0.01076609 0.003853774
#                    robustCV
# S000004-37958-003 0.2956749
# S000004-38071-002 0.2886340
# S000004-38094-003 0.2998139
# S000004-38120-002 0.3165325
# S000004-38134-002 0.3129712
# S000004-38138-003 0.3579549
```

## Sample Quality Index

Calculate a mean coverage depth and a window coefficient of variation

``` r
ptm=proc.time()[3]
MCD.mat = get_MCD(genelist, pileupPath, sampleInfo)
proc.time()[3]-ptm
# elapsed 
#  17.593

ptm=proc.time()[3]
wCV.mat = get_wCV(genelist, pileupPath, sampleInfo, rnum=100, method=1, winSize=20, egPct=10)
proc.time()[3]-ptm
# elapsed 
#  78.392 
```

Sample quality index

``` r
ptm=proc.time()[3]
result = get_SQI(MCD=MCD.mat, wCV=wCV.mat, rstPct=20, obsPct=50)
proc.time()[3]-ptm
# elapsed 
#   2.512 

auc.vec <- result$auc.vec
table(auc.vec$SQI)
 # Bad Good 
 #  19  152 

plot_SQI(SQIresult=result)
```
![](figures/Allianceex_SQI_v058.png)<!-- -->
<div align="center">
  <img width="70%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/figures/Allianceex_SQI_v058.png">
</div>

### Update gene body coverage with good quality samples

``` r
GBCg0 = plot_GBCg(stat=2, plot=TRUE, sampleInfo, GBCresult=GBC0, auc.vec=result$auc.vec)
GBCg5 = plot_GBCg(stat=2, plot=TRUE, sampleInfo, GBCresult=GBC5, auc.vec=result$auc.vec)

pg0 <- GBCg0$plotRI +
  coord_cartesian(ylim=c(0, 0.017)) +
  ggtitle("0~5 kb")

pg5 <- GBCg5$plotRI +
  coord_cartesian(ylim=c(0, 0.017)) +
  ggtitle("5+ kb")

ggpubr::ggarrange(pg0, pg5, common.legend=TRUE, legend="bottom", nrow=1)


pg0 <- GBCg0$plotPD +
  coord_cartesian(ylim=c(0, 0.017)) +
  ggtitle("0~5 kb")

pg5 <- GBCg5$plotPD +
  coord_cartesian(ylim=c(0, 0.017)) +
  ggtitle("5+ kb")

ggpubr::ggarrange(pg0, pg5, common.legend=TRUE, legend="bottom", nrow=1)
```
![](figures/Allianceex_GBCgRI_v058.png)<!-- -->
<div align="center">
  <img width="70%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/figures/Allianceex_GBCgRI_v058.png">
</div>

![](figures/Allianceex_GBCgPD_v058.png)<!-- -->
<div align="center">
  <img width="70%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/figures/Allianceex_GBCgPD_v058.png">
</div>

## Principal Component Analysis

``` r
# Convert NA to 0
wCV.mat2 <- wCV.mat
wCV.mat2[is.na(wCV.mat2)] <- 0

# Remove constant rows (zero variance rows)
wCV.mat2 <- wCV.mat2[apply(wCV.mat2, 1, var)!=0, ]
dim(wCV.mat2) # 1000  171

# Perform PCA
pca_result <- prcomp(t(wCV.mat2), scale=TRUE)

# Contributions of each gene to the PC1
pc1_contributions <- abs(pca_result$rotation[, 1])
top_genes <- order(pc1_contributions, decreasing=TRUE)
```

## window CV Heatmap
![](figures/Allianceex_wCVheatmap_v058-4.png)<!-- -->
<div align="center">
  <img width="75%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/figures/Allianceex_wCVheatmap_v058-4.png">
</div>
