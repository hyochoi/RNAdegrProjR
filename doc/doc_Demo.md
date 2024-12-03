Demo
================

This Alliance example consists of 101 selected genes and FFT 171
samples. The union transcript is used to extract only exon pileup.

## Gene Length Normalization

Build coverage pileup from raw pileup (part_intron) to pileupData
(only_exon)

``` r
get_pileupData <- function(g) {
  load(file=rawPath[g])
  pileupData = SCISSOR::build_pileup(Pileup=pileup, regions=regions, inputType="part_intron", outputType="only_exon")
  colnames(pileupData) <- colnames(pileup) # to keep the original sample IDs

  save(pileupData, file=paste0(folder_newpath, sprintf("%s_pileupData_only_exon.RData", genes[g])))
}

# lapply(1:length(genes), get_pileupData)
```

Normalization by methods and results from the g-th gene of i-th sample

``` r
g=50; i=100
genes[g]; sampleInfo$SampleID[i]
```

    ## [1] "ERO1A"

    ## [1] "S000004-38105-003"

``` r
load(file=rawPath[g])
pileupData = SCISSOR::build_pileup(Pileup=pileup, regions=regions, inputType="part_intron", outputType="only_exon")
geneRanges = get_Ranges(Gene=Ranges$Gene, regions=regions, outputType="only_exon")
colnames(pileupData) <- colnames(pileup) # to keep the original sample IDs

# method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively.
norm1 = norm_pileup.gene(pileupData, rnum=100, method=1)
norm2 = norm_pileup.gene(pileupData, rnum=100, method=2)


par(mfrow=c(1,3))

# Original plot_pileup
SCISSOR::plot_pileup(Pileup=log10(pileupData+1), case=i, Ranges=geneRanges, main=sampleInfo$SampleID[i])

# Methods
plot(log10(norm1[,i]+1), type='l', main="Method 1: Raw value", ylab='')
plot(log10(norm2[,i]+1), type='l', main="Method 2: Interpolation", ylab='')
```

![](figures/unnamed-chunk-5-1.png)<!-- -->

``` r
list1 = norm_pileup.list(pileupPath, geneNames=genes, rnum=100, method=1)
list2 = norm_pileup.list(pileupPath, geneNames=genes, rnum=100, method=2)
cbind(list1[[g]][91:100, i], list2[[g]][91:100, i])
```

    ##     [,1]        [,2]
    ## 91   733  825.508318
    ## 92   722  733.766630
    ## 93  1099  887.586518
    ## 94  1083 1064.850365
    ## 95   439  477.522727
    ## 96   109  169.293864
    ## 97   150  155.923548
    ## 98   190  157.173955
    ## 99    99   57.557664
    ## 100    5    6.348469

Metrics from scaled normalized transcript coverage by margin

``` r
# margin 1, 2, and 3 return metrics per sample, per gene, and across the genes per sample, respectively.
met1 = get_metrics(pileupPath, geneNames=genes, rnum=100, method=1, scale=TRUE, margin=1)
met2 = get_metrics(pileupPath, geneNames=genes, rnum=100, method=1, scale=TRUE, margin=2)
met3 = get_metrics(pileupPath, geneNames=genes, rnum=100, method=1, scale=TRUE, margin=3)
head(met1)
```

    ##                         mean          sd        CV     median         mad
    ## S000004-37958-003 0.01071488 0.008166250 0.7621408 0.01067084 0.004007677
    ## S000004-38071-002 0.01059379 0.005078067 0.4793437 0.01062222 0.003744951
    ## S000004-38094-003 0.01054652 0.010619150 1.0068865 0.01050303 0.003921993
    ## S000004-38120-002 0.01060603 0.004175746 0.3937144 0.01066341 0.003861760
    ## S000004-38134-002 0.01046946 0.003971038 0.3792971 0.01056412 0.003678772
    ## S000004-38138-003 0.01079762 0.005678570 0.5259095 0.01069435 0.004303059
    ##                    robustCV
    ## S000004-37958-003 0.3755728
    ## S000004-38071-002 0.3525582
    ## S000004-38094-003 0.3734155
    ## S000004-38120-002 0.3621508
    ## S000004-38134-002 0.3482327
    ## S000004-38138-003 0.4023676

``` r
head(met2)
```

    ##             mean          sd        CV      median         mad  robustCV
    ## A2M   0.01048835 0.003827195 0.3648997 0.011934423 0.002275591 0.1906746
    ## ABCC2 0.01519139 0.008953634 0.5893886 0.013398457 0.005672699 0.4233845
    ## ACTB  0.01208595 0.006605041 0.5465056 0.012155400 0.009214851 0.7580870
    ## ACTG1 0.01026742 0.004109046 0.4002022 0.012004023 0.003147533 0.2622065
    ## AHNAK 0.01009524 0.002201511 0.2180740 0.010613081 0.001621049 0.1527406
    ## ANXA1 0.01035861 0.004519427 0.4362966 0.009367601 0.003959851 0.4227178

``` r
head(met3)
```

    ##              sample gene        mean          sd        CV     median
    ## 1 S000004-37958-003  A2M 0.009999675 0.003109082 0.3109183 0.01136437
    ## 2 S000004-38071-002  A2M 0.009999673 0.002882120 0.2882214 0.01120771
    ## 3 S000004-38094-003  A2M 0.009999635 0.002932358 0.2932465 0.01118402
    ## 4 S000004-38120-002  A2M 0.009999687 0.003253696 0.3253798 0.01141832
    ## 5 S000004-38134-002  A2M 0.009999657 0.003364398 0.3364513 0.01166245
    ## 6 S000004-38138-003  A2M 0.009999660 0.003542085 0.3542205 0.01151411
    ##           mad  robustCV
    ## 1 0.001697565 0.1493760
    ## 2 0.001841086 0.1642696
    ## 3 0.001847067 0.1651524
    ## 4 0.001801976 0.1578145
    ## 5 0.001662940 0.1425892
    ## 6 0.002039569 0.1771364

# Sample Quality Index (SQI)

Calculate a mean coverage depth and a window coefficient of variation

``` r
MCD.mat = get_MCD(genelist, pileupPath, sampleInfo)
wCV.mat = get_wCV(genelist, pileupPath, sampleInfo, rnum=100, method=1, winSize=20, egPct=10)
```

Sample quality index

``` r
result = get_SQI(MCD=MCD.mat, wCV=wCV.mat, rstPct=10, logMCD=TRUE)
head(result$auc.vec)
```

    ## # A tibble: 6 × 6
    ##   Sample              AUC  zScore      PD SQI.zScore SQI.PD
    ##   <chr>             <dbl>   <dbl>   <dbl> <chr>      <chr> 
    ## 1 S000004-37933-003  2.50  1.86    1.66   Good       Good  
    ## 2 S000004-37935-002  2.21  0.786   0.701  Good       Good  
    ## 3 S000004-37937-002  2.47  1.74    1.55   Good       Good  
    ## 4 S000004-37939-002  2.40  1.49    1.33   Good       Good  
    ## 5 S000004-37941-002  1.97 -0.0857 -0.0948 Good       Good  
    ## 6 S000004-37943-002  2.23  0.862   0.769  Good       Good

Fitted wCV curves and quality detection using robust z-score and
projection depth

![](figures/unnamed-chunk-14-1.png)<!-- -->
