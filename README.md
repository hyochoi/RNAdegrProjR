# RNAdegrProjR
Deciphering RNA degradation: Insights from a comparative analysis of paired fresh frozen/FFPE total RNA-seq


## RNAdegr Description
`RNAdegr` investigates RNA degradation patterns in fresh frozen mRNA-seq [FFM], fresh frozen total RNA-seq [FFT], and FFPE total RNA-seq [PET]. `RNAdegr` measures noise patterns in selected samples by a method called windowCV (wCV) utilizing the coefficient of variance (CV) along the transcript length. The simple measure dividing the expression value of other proteins with certain lncRNA or mtRNA inferring the degree of RNA degradation is also included.


## Data Pre-processing
- Gene information table
- Sample information table
- Gene expression normalization[^1]: FPKM and TPM[^2], FPKM-UQ[^3], TMM[^4]
[^1]: https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html
[^2]: https://doi.org/10.1186/s12967-021-02936-w
[^3]: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
[^4]: https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html


## Documentation
- R function overview
- TCGA example


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
> [!NOTE]
> `RNAdegr` requires the `SCISSOR` package, and its dependencies are at https://github.com/hyochoi/SCISSOR.


## References
