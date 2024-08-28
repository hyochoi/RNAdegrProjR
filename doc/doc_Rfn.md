# R Function Overview


## Coverage Pileup


## Gene Length Normalization
- `norm_pileup.spl`: gene length normalization for a single pileup (for gene 1, sample 1)
- `norm_pileup.gene`: gene length normalization for a pileup matrix (for gene 1, all samples)
- `norm_pileup.list`: gene length normalization for pileup lists (for all genes, all samples)
- `scale_pileup.list`: scale normalized transcript coverage using gene length normalization


### Method options
- Method 1: find read depth at even points (green points)
- Method 2: find geometric mean (blue points) using read depth at odd points (red points)
![alt text](https://github.com/hyochoi/RNAdegrProjR/blob/main/images/norm_pileup_methods2.png?raw=true)
  - We recommend using the gene length normalization methods for gene length is at least 201 (2*the number of regions+1) where the number of regions is 100.


## windowCV
- `get_windowCV`
