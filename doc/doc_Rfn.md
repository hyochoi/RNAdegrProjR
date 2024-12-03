# R Functions


## Gene Length Normalization
- `norm_pileup.spl`: gene length normalization for a single pileup (for gene 1, sample 1).
- `norm_pileup.gene`: gene length normalization for a pileup matrix (for gene 1, all samples).
- `norm_pileup.list`: gene length normalization for pileup lists (for all genes, all samples).
- `scale_pileup.list`: scale normalized transcript coverage using gene length normalization.
- `get_metrics`: get metrics from scaled normalized transcript coverage.

#### Method options
- Method 1: find read depth at even points (green points).
- Method 2: find geometric mean (blue points) using read depth at odd points (red points).
  - We recommend using the gene length normalization methods for gene length is at least 201 (2*the number of regions+1) where the number of regions is 100.
<div align="center">
  <img width="60%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/images/norm_pileup_methods2.png">
</div>


## Sample Quality Index
- `get_focusPileup`: get a focused pileup of exon location for the g-th gene.
- `get_MCD`: get a mean coverage depth (MCD) for genes and samples.
- `get_wCV`: get a window coefficient of variation (wCV) for genes and samples.
- `get_SQI`: get a sample quality index (SQI) for samples.
  
