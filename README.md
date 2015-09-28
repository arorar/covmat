# covmat
Package Development for GSOC 2015

[![Build Status](https://travis-ci.org/arorar/covmat.svg?branch=master)](https://travis-ci.org/arorar/covmat) 
[![Coverage Status](https://coveralls.io/repos/arorar/covmat/badge.svg?branch=master)](https://coveralls.io/github/arorar/covmat)
[![CRAN](http://www.r-pkg.org/badges/version/covmat)](http://cran.rstudio.com/package=covmat) [![Downloads](http://cranlogs.r-pkg.org/badges/covmat?color=brightgreen)](http://www.r-pkg.org/pkg/covmat)

Covmat is a collection of techniques for estimating convariance matrices. Covariance matrices can be built using missing data. Stambaugh Estimation and FMMC methods can be used to construct such matrices. Covariance matrices can be built by denoising or shrinking the eigenvalues of a sample covariance matrix. Such techniques work by exploiting the tools in Random Matrix Theory to analyse the distribution of eigenvalues. Covariance matrices can also be built assuming that data has many underlying regimes. Each regime is allowed to follow a Dynamic Conditional Correlation model. Robust covariance matrices can be constructed by multivariate cleaning and smoothing of noisy data.

Installation
------------

To get started, you can install the package from github using `devtools`.

``` r
library(devtools)
install_github("arorar/covmat")
```

Examples
--------

Detailed information on covmat's functionality and use can be found by reading the **[covmat vignette](https://github.com/arorar/covmat/blob/master/inst/doc/CovarianceEstimation.pdf)**
