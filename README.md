
# scAMACE (R implementation)

scAMACE (integrative Analysis of single-cell Methylation, chromatin ACcessibility, and gene Expression)

A model-based approach to the joint analysis of single-cell data on chromatin accessibility, gene expression and methylation.

## 1. Installation

You can install the released version of scAMACE from Github:

``` r
library(devtools)
devtools::install_github("cuhklinlab/scAMACE")
```

Package 'betareg' is also required for the implementation of beta regression:

```{r}
install.packages('betareg')
```


## 2. Main Functions

`cal_M_step`: scAMACE expectation-maximization (EM) implementation, a model-based approach to the joint clustering of single-cell data on chromatin accessibility, gene expression and methylation.

`cal_post`: Calculate the posterior probability for one iteration in the EM algorithm.

`cal_E_rna`: Perform E-step (i.e. calculate the expectations of missing data) for one iteration of scRNA-Seq or sc-methylation data in the EM algorithm.

`cal_E_acc`: Perform E-step (i.e. calculate the expectations of missing data) for one iteration of scCAS data in the EM algorithm.

`simData_3data`: Generate simulation data x, y and z.

## 3. Datasets and Examples

Please refer to the [vigenette](https://github.com/cuhklinlab/scAMACE/blob/main/vignette/vignette.md) with several examples for a quick guide to scAMACE package.

## 4. Reference
Jiaxuan Wangwu, Zexuan Sun, Zhixiang Lin: [scAMACE: Model-based approach to the joint analysis of single-cell data on chromatin accessibility, gene expression and methylation.](https://www.biorxiv.org/content/10.1101/2021.03.29.437485v1)

