# RNAdegrProjR

Despite easy storage and cost-effectiveness advantages, formalin-fixed paraffin-embedded (FFPE) samples have the disadvantage of inevitable chemical-induced RNA degradation. While the 3' bias due to the characteristics of the mRNA-seq platform allows the measure of RNA degradation levels in mRNA-seq data, there is still no clear measure for RNA degradation for total RNA-seq and FFPE samples. `RNAdegrProjR` investigates RNA degradation patterns in fresh frozen mRNA-seq [FFM], fresh frozen total RNA-seq [FFT], and FFPE total RNA-seq [PET], which measures noise patterns in total RNA-seq and FFPE samples by a method called windowCV (wCV). The sample quality index can be achieved by outlier detection from the relation between wCV and mean coverage depth.

- Abstract: https://doi.org/10.1158/1538-7445.AM2024-2323


## Documentation

- [R functions](https://github.com/hyochoi/RNAdegrProjR/blob/master/doc/doc_Rfn.md)
- [Demo](https://github.com/hyochoi/RNAdegrProjR/blob/master/doc/doc_Demo.html)


## Data Processing

- Gene information table: combine exon, intron, and transcript level datasets to gene level data.
- Sample information table: summarize RNA-seq platform, metrics about the alignment of RNA-seq reads, and degradation indexes.
- Gene expression transformation: normalize raw read counts using FPKM, FPKM-UQ, TPM, and TMM methods
<!---[^1]: https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html --->
<!---[^2]: https://doi.org/10.1186/s12967-021-02936-w --->
<!---[^3]: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/ --->
<!---[^4]: https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html --->


## Data Analysis

- FFT vs. PET mean scatter and violin plots to compare lncRNA and protein-coding
- Pileup plots to investigate the difference in read depth by sample and RNA-seq platform
- Genome alignment profiles using metrics about the alignment of RNA-seq reads
- Normalized transcript coverage by degradation indexes
- Gene length normalization steps for gene KEAP1
- GC percentage TMM scatter plot
- Ratio heatmap with hierarchical clustering of genes and patients
<!---[^5]: https://doi.org/10.1186/s12864-017-3827-y --->
<!---[^6]: https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics --->
<!---[^7]: https://doi.org/10.1186/1471-2164-15-419 --->
