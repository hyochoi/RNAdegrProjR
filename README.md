# RNAdegrProjR

Despite easy storage and cost-effectiveness advantages, formalin-fixed paraffin-embedded (FFPE) samples have the disadvantage of inevitable chemical-induced RNA degradation. While the 3' bias due to the characteristics of the mRNA-seq platform allows the measure of RNA degradation levels in mRNA-seq data, there is still no clear measure for RNA degradation for total RNA-seq and FFPE samples. `RNAdegrProjR` investigates RNA degradation patterns in fresh frozen mRNA-seq [FFM], fresh frozen total RNA-seq [FFT], and FFPE total RNA-seq [PET], which measures noise patterns in total RNA-seq and FFPE samples by a method called windowCV (wCV). The sample quality index can be achieved by outlier detection from the relation between wCV and mean coverage depth.

- Abstract: https://doi.org/10.1158/1538-7445.AM2024-2323

<div align="center">
  <img width="75%" src="https://github.com/hyochoi/RNAdegrProjR/blob/main/figures/Overview.png">
</div>


## Documentation
The documentation and demo are available at [RNAdegrDocs](https://bookdown.org/sqr_yeon/RNAdegrDocs/).


## Citation
Choi, W., Yeon, M., Lee, J., Choi, H., Hayes, D.N.(2025+).


## Installation
- From GitHub
  ```r
  if (!requireNamespace("devtools", quietly=TRUE)) {install.packages("devtools")}
  devtools::install_github("hyochoi/RNAdegrProjR", dependencies=TRUE)
  ```

- From Bioconductor
  ```r
  if (!requireNamespace("BiocManager", quietly=TRUE)) {install.packages("BiocManager")}
  BiocManager::install("RNAdegrProjR", dependencies=TRUE)
  ```
