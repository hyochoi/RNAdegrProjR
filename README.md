# RNAdegrProjR
Comparison and assessment of mRNA degradation in fresh frozen tissue and formalin-fixed paraffin-embedded tissue


## Abstract
In this study, we examine RNA-seq data paired in fresh frozen (FF) and FFPE from multiple tumor types obtained from a subset of The Cancer Genome Atlas (TCGA) cohorts. We compare FF mRNA-seq (FFM), FF total RNA-seq (FFT), and FFPE total RNA-seq (PET) and evaluate two different RNA preservation methods (FF vs. FFPE) as well as two different RNA preparation protocols (mRNA-seq vs. Ribo-Zero-seq). First, we perform genome-wide comparisons of gene expression profiles with all three pairs of FFM, FFT, and PET using unsupervised/supervised clustering methods. Second, we evaluate gene-specific and sample-specific degradation using SCISSOR and identify severely degraded samples. Our results demonstrate that the impact of formalin-induced degradation is highly gene-specific that might be associated with multiple factors such as RNA localization, gene length, and GC concentration. We also show that the level of formalin-induced degradation could be different in FFPE samples and thus should be correctly adjusted for the downstream gene expression analysis.


## Data pre-processing
1. Gene information table
1. Sample information table
1. Gene expression normalization[^1]: FPKM and TPM[^2], FPKM-UQ[^3], TMM[^4]
[^1]: https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html
[^2]: https://doi.org/10.1186/s12967-021-02936-w
[^3]: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
[^4]: https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html


## R functions
### Gene length normalization
- Method 1: find read depth at even points (green points)
- Method 2: find geometric mean (blue points) using read depth at odd points (red points)
![alt text](https://github.com/hyochoi/RNAdegrProjR/blob/main/images/norm_pileup_methods2.png?raw=true)
> [!NOTE]
> We recommend using the gene length normalization methods for gene length is at least 201 (2*the number of regions+1) where the number of regions is 100.


## Analysis
- Pileup plot for each sample (The difference in read depth between samples)
- Normalized transcript coverage[^5] using Picard tool CollectRnaSeqMetrics[^6]
- Genome alignment profiles[^7]
- FFT vs. PET mean scatter plot (Comparison between lncRNA and protein-coding)
- Ratio heatmap with ward.D clustering of genes and patients
[^5]: https://doi.org/10.1186/s12864-017-3827-y
[^6]: https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics
[^7]: https://doi.org/10.1186/1471-2164-15-419


## Installation
```
if ("devtools" %in% rownames(installed.packages()) == FALSE) {install.packages("devtools")}
library(devtools)

install_github("hyochoi/RNAdegr")
library(RNAdegr)
```


## References
