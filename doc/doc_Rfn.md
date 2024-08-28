# R Function Overview


## Step 0. Coverage Pileup


## Step 1. Gene Length Normalization
- `norm_pileup.spl`: gene length normalization for a single pileup (for gene 1, sample 1)
- `norm_pileup.gene`: gene length normalization for a pileup matrix (for gene 1, all samples)
- `norm_pileup.list`: gene length normalization for pileup lists (for all genes, all samples)
- `scale_pileup.list`: scale normalized transcript coverage using gene length normalization

### Method options
- Method 1: find read depth at even points (green points)
- Method 2: find geometric mean (blue points) using read depth at odd points (red points)
  - We recommend using the gene length normalization methods for gene length is at least 201 (2*the number of regions+1) where the number of regions is 100.
<div align="center">
  <img width="70%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/images/norm_pileup_methods2.png">
</div>


## Step 2. Gene Expression Normalization[^1]
- FPKM and TPM[^2], FPKM-UQ[^3], TMM[^4]
[^1]: https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html
[^2]: https://doi.org/10.1186/s12967-021-02936-w
[^3]: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
[^4]: https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html


## Step 3. windowCV
- `get_windowCV`
