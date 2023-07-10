## -------------------------------------------------
##    code to prepare `DATASET` dataset goes here
## -------------------------------------------------

## ---- output


## ---- renv

root <- "."
devtools::load_all(root, quiet = TRUE)

## ---- prepare the workspace

# imports
suppressPackageStartupMessages({
  library(dplyr) # For manipulating the data
  library(data.table) # For working with data.tables
  library(biomaRt) # For getting gene information
  # library(org.Hs.eg.db) # For extracting entrezIDs
  library(edgeR) # For TMM normalization
})

## ---- 1. Get gene information table

# Miyeon - Could you add your code to create the variables of GeneInfo.csv from the raw datasets?

geneinfo.csv <- read.csv("data-raw/GeneInfo.csv")
geneinfo.all <- geneinfo.csv[,-1]

## ---- 2. Get sample information table



## ---- 3. Gene expression normalization

# Miyeon - Could you add your code for normalization?



## --------------------------------------
##     save processed data to data/
## --------------------------------------

usethis::use_data(geneinfo.all, overwrite = TRUE)   # gene information table for all

