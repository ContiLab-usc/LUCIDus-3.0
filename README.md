
# LUCIDus: Integreted clustering with multi-view data Version 3.0.1

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/LUCIDus?color=green)](https://cran.r-project.org/package=LUCIDus)
![](https://cranlogs.r-pkg.org/badges/grand-total/LUCIDus?color=blue)
[![](https://raw.githubusercontent.com/USCbiostats/badges/master/tommy-image-badge.svg)](https://image.usc.edu)
<!-- badges: end -->



The **LUCIDus** package implements the statistical method LUCID proposed
in the research paper [A Latent Unknown Clustering Integrating
Multi-Omics Data (LUCID) with Phenotypic
Traits](https://doi.org/10.1093/bioinformatics/btz667)
(*Bioinformatics*, 2020). LUCID conducts integrated clustering by using
multi-view data, including exposures, omics data with/without outcome.
**LUCIDus** features variable selection, incorporating missingness in
omics data, visualization of LUCID model via Sankey diagram, bootstrap
inference and functions for tuning model parameters.

If you are interested in integration of omic data to estiamte mediator
or latent structures, please check out [Conti
Lab](https://contilab.usc.edu/about/) to learn more.

LUCID version 3, a major update and enhancement from the original release. LUCID version 3 is more robust, and features more powerful model selection, model visualization, inference based on bootstrap resampling, and implements different analysis strategies for multi-omics data with multiple layers. It also incorporates methods to deal with missingness in multi-omics data. 
![plot](./figure/fig1.png)

## Installation

You can install the development version of LUCIDusM from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ContiLab-usc/LUCIDus-3.0",ref="main",auth_token = "xxx")
```
Note that this repo is now private, so only authorized users can download this package. Please go to [tokens](https://github.com/settings/tokens) to obtain your personal authorized token and input it into auth_token = "xxx" to download this package.

## Workflow
![plot](./figure/fig2.png)
## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(LUCIDus)
## basic example code
```

