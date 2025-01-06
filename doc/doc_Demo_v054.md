Demo
================

## Data Description

This Alliance example consists of 1,000 selected genes and FFT 171
samples. The union transcript is used to extract only exon pileup. We
need datasets such as `genelist`, coverage `pileup` (located in
`pileupPath`), and `sampleInfo` to obtain the sample quality index
outputs and plots. `geneInfo` is optional in case you want to compare
results by properties of genes.

To keep only exon location, we first build coverage `pileup` from raw
pileup (part_intron) to `pileupData` (only_exon). Let’s take a look into
the process of the first and the last genes. The equal number of
positions will be selected from the different genomic positions:
`LINC01772` and `MIR133A1HG` have 3,245 and 5,825 positions,
respectively.

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
![](figures/Allianceex_GAP_v054.png)<!-- -->
<div align="center">
  <img width="90%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/figures/Allianceex_GAP_v054.png">
</div>

## Gene Body Coverage

See R function in Documentation for details. Mean scaled normalized
coverage at the gene body percentile for A. all genes; B. 0~5 kb; C. 5+
kb

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


GBC = plot_GBC(pileupPath2, geneNames=genelist2, rnum=100, method=1, scale=TRUE, stat=2, plot=TRUE, sampleInfo)
GBC0 = plot_GBC(pileupPath2Len0, geneNames=genelist2Len0, rnum=100, method=1, scale=TRUE, stat=2, plot=TRUE, sampleInfo)
GBC5 = plot_GBC(pileupPath2Len5, geneNames=genelist2Len5, rnum=100, method=1, scale=TRUE, stat=2, plot=TRUE, sampleInfo)

p <- GBC$plot +
  coord_cartesian(ylim=c(0, 0.017))
p0 <- GBC0$plot +
  coord_cartesian(ylim=c(0, 0.017))
p5 <- GBC5$plot +
  coord_cartesian(ylim=c(0, 0.017))

ggpubr::ggarrange(p, p0, p5, labels=c("A", "B", "C"), common.legend=TRUE, legend="bottom", nrow=1)
```
![](figures/Allianceex_GBC_v054.png)<!-- -->
<div align="center">
  <img width="90%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/figures/Allianceex_GBC_v054.png">
</div>

``` r
## Metrics from scaled normalized transcript coverage by margin
# margin 1, 2, and 3 return metrics per sample, per gene, and across the genes per sample, respectively.
# met1 = get_metrics(pileupPath, geneNames=genes, rnum=100, method=1, scale=TRUE, margin=1)
# met2 = get_metrics(pileupPath, geneNames=genes, rnum=100, method=1, scale=TRUE, margin=2)
# met3 = get_metrics(pileupPath, geneNames=genes, rnum=100, method=1, scale=TRUE, margin=3)
# head(met1)
# head(met2)
# head(met3)
```

## Sample Quality Index (SQI)

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
![](figures/Allianceex_SQI_v054.png)<!-- -->
<div align="center">
  <img width="90%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/figures/Allianceex_SQI_v054.png">
</div>

Update gene body coverage with good quality samples

``` r
sampleInfo2 <- sampleInfo[match(c(auc.vec[auc.vec$SQI=="Good", c("Sample")])$Sample, sampleInfo$SampleID), ]


pg <- plot_GBCg(stat=2, plot=TRUE, sampleInfo, GBCresult=GBC, auc.vec=result$auc.vec)$plot +
  coord_cartesian(ylim=c(0, 0.017))
p0g <- plot_GBCg(stat=2, plot=TRUE, sampleInfo, GBCresult=GBC0, auc.vec=result$auc.vec)$plot +
  coord_cartesian(ylim=c(0, 0.017))
p5g <- plot_GBCg(stat=2, plot=TRUE, sampleInfo, GBCresult=GBC5, auc.vec=result$auc.vec)$plot +
  coord_cartesian(ylim=c(0, 0.017))

ggpubr::ggarrange(pg, p0g, p5g, labels=c("D", "E", "F"), common.legend=TRUE, legend="bottom", nrow=1)
```
![](figures/Allianceex_GBCg_v054.png)<!-- -->
<div align="center">
  <img width="90%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/figures/Allianceex_GBCg_v054.png">
</div>
