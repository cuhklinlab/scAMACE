
# scAMACE (R implementation)

scAMACE (integrative Analysis of single-cell Methylation, chromatin ACcessibility, and gene Expression)

A model-based approach to the joint analysis of single-cell data on chromatin accessibility, gene expression and methylation.

## 1. Installation

You can install the released version of scAMACE from Github:

``` r
library(devtools)
devtools::install_github("WWJiaxuan/scAMACE")
```

## 2. Main Functions

`cal_M_step`: scAMACE expectation-maximization (EM) implementation, a model-based approach to the joint clustering of single-cell data on chromatin accessibility, gene expression and methylation.

`cal_post`: Calculate the posterior probability for one iteration in the EM algorithm.

`cal_E_rna`: Perform E-step (i.e. calculate the expectations of missing data) for one iteration of scRNA-Seq or sc-methylation data in the EM algorithm.

`cal_E_acc`: Perform E-step (i.e. calculate the expectations of missing data) for one iteration of scCAS data in the EM algorithm.

## 3. Datasets and Examples

Please refer to the [vigenette](https://github.com/WWJiaxuan/scAMACE/tree/main/vignette) with several examples for a quick guide to scAMACE package.
