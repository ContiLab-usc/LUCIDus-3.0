
# LUCIDus Version 3.0.1

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/LUCIDus?color=green)](https://cran.r-project.org/package=LUCIDus)
![](https://cranlogs.r-pkg.org/badges/grand-total/LUCIDus?color=blue)
[![](https://raw.githubusercontent.com/USCbiostats/badges/master/tommy-image-badge.svg)](https://image.usc.edu)
<!-- badges: end -->

![plot](./figure/fig1.png)

LUCID version 3, a major update and enhancement from the original release. LUCID version 3 is more robust, and features more powerful model selection, model visualization, inference based on bootstrap resampling, and implements different analysis strategies for multi-omics data with multiple layers. It also incorporates methods to deal with missingness in multi-omics data. 
![plot](./figure/fig2.png)

## Installation

You can install the development version of LUCIDusM from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ContiLab-usc/LUCIDus-3.0",ref="main",auth_token = "xxx")
```
Note that this repo is now private, so only authorized users can download this package. Please go to [tokens](https://github.com/settings/tokens) to obtain your personal authorized token and input it into auth_token = "xxx" to download this package.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(LUCIDus)
## basic example code
```

