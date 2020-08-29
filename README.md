
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ripserr: Calculate Persistent Homology of Vietoris-Rips and Cubical Complexes using Ripser in R

[![Travis-CI Build
Status](https://travis-ci.org/rrrlw/ripserr.svg?branch=master)](https://travis-ci.org/rrrlw/ripserr)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/rrrlw/ripserr?branch=master&svg=true)](https://ci.appveyor.com/project/rrrlw/ripserr)
[![Codecov test
coverage](https://codecov.io/gh/rrrlw/ripserr/branch/master/graph/badge.svg)](https://codecov.io/gh/rrrlw/ripserr?branch=master)

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN
version](http://www.r-pkg.org/badges/version/ripserr)](https://CRAN.R-project.org/package=ripserr)
[![CRAN
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/ripserr)](https://CRAN.R-project.org/package=ripserr)

## Overview

ripserr ports the [Ripser](https://arxiv.org/1908.02518) and [Cubical
Ripser](https://arxiv.org/abs/2005.12692) persistent homology
calculation engines from C++ via
[Rcpp](https://CRAN.R-project.org/package=Rcpp). It can be used as a
convenient and rapid calculation tool in topological data analysis
pipelines.

## Installation

``` r
# install development version
devtools::install_github("rrrlw/ripserr")

# install from CRAN (not available yet)
# install.packages("ripserr")
```

## Sample code

Ripser (Vietoris-Rips complex) can be used as follows for data with
dimension greater than or equal to 2.

``` r
# load ripserr
library("ripserr")

set.seed(42)
SIZE <- 100

# 2-dimensional example
dataset2 <- rnorm(SIZE * 2)
dim(dataset2) <- c(SIZE, 2)
vr_phom2 <- vietoris_rips(dataset2)
head(vr_phom2)
#>   dimension birth      death
#> 1         0     0 0.01004861
#> 2         0     0 0.02923702
#> 3         0     0 0.04550504
#> 4         0     0 0.06829826
#> 5         0     0 0.06853393
#> 6         0     0 0.07187663
tail(vr_phom2)
#>     dimension     birth     death
#> 112         1 0.3916344 0.4239412
#> 113         1 0.3906770 0.5577989
#> 114         1 0.3880186 0.4029842
#> 115         1 0.3703398 0.5007012
#> 116         1 0.3330234 0.3416054
#> 117         1 0.2418318 0.2504820

# 3-dimensional example
dataset3 <- rnorm(SIZE * 3)
dim(dataset3) <- c(SIZE, 3)
vr_phom3 <- vietoris_rips(dataset3, dim = 2) # default: dim = 1
head(vr_phom3)
#>   dimension birth     death
#> 1         0     0 0.1282935
#> 2         0     0 0.1421812
#> 3         0     0 0.1516424
#> 4         0     0 0.1819928
#> 5         0     0 0.1858051
#> 6         0     0 0.2114116
tail(vr_phom3)
#>     dimension     birth     death
#> 132         1 0.5212961 0.5233529
#> 133         2 1.1829207 1.1999911
#> 134         2 1.1194325 1.3245908
#> 135         2 1.0707410 1.0914850
#> 136         2 0.9433034 0.9867254
#> 137         2 0.6882204 0.6913078
```

Cubical Ripser (cubical complex) can be used as follows for data with
dimension equal to 2, 3, or 4.

``` r
# load ripserr
library("ripserr")

set.seed(42)
SIZE <- 10

# 2-dimensional example
dataset2 <- rnorm(SIZE ^ 2)
dim(dataset2) <- rep(SIZE, 2)
cub_phom2 <- cubical(dataset2)
head(cub_phom2)
#>   dimension      birth      death
#> 1         0 -1.1943289 -0.8607926
#> 2         0 -2.4142076 -0.8509076
#> 3         0 -0.8113932 -0.7844590
#> 4         0 -1.7170087 -0.7844590
#> 5         0 -0.7272921 -0.5428288
#> 6         0 -0.9535234 -0.5428288
tail(cub_phom2)
#>    dimension     birth     death
#> 23         1 0.8217731 0.9333463
#> 24         1 0.7681787 1.0385061
#> 25         1 0.7581632 1.5757275
#> 26         1 0.7208782 1.3025426
#> 27         1 0.6792888 1.4441013
#> 28         1 0.6359504 1.8951935

# 3-dimensional example
dataset3 <- rnorm(SIZE ^ 3)
dim(dataset3) <- rep(SIZE, 3)
cub_phom3 <- cubical(dataset3)
head(cub_phom3)
#>   dimension     birth     death
#> 1         0 -1.926167 -1.737728
#> 2         0 -1.737297 -1.439229
#> 3         0 -1.924950 -1.439229
#> 4         0 -1.500221 -1.354600
#> 5         0 -2.277778 -1.354600
#> 6         0 -1.682481 -1.306676
tail(cub_phom3)
#>     dimension     birth    death
#> 325         2 1.2488637 1.258482
#> 326         2 1.2009654 2.036972
#> 327         2 1.0452759 1.199978
#> 328         2 0.9885968 1.809382
#> 329         2 0.9310749 1.179696
#> 330         2 0.8447922 1.709689

# 4-dimensional example
dataset4 <- rnorm(SIZE ^ 4)
dim(dataset4) <- rep(SIZE, 4)
cub_phom4 <- cubical(dataset4)
head(cub_phom4)
#>   dimension     birth     death
#> 1         0 -1.986299 -1.923519
#> 2         0 -1.822606 -1.816506
#> 3         0 -1.776392 -1.710786
#> 4         0 -1.833663 -1.710387
#> 5         0 -1.947054 -1.704791
#> 6         0 -1.701462 -1.639160
tail(cub_phom4)
#>      dimension    birth    death
#> 4330         3 1.676609 2.019277
#> 4331         3 1.675766 1.932152
#> 4332         3 1.669449 2.149646
#> 4333         3 1.662486 1.863734
#> 4334         3 1.535361 1.963609
#> 4335         3 1.349235 2.263581
```

## Functionality

1.  Calculation of persistent homology of Vietoris-Rips complexes using
    Ripser (function named `vietoris_rips`).
2.  Calculation of persistent homology of cubical complexes using
    Cubical Ripser (function named `cubical`).

## Citation

If you use the ripserr package in your work, please consider citing the
following (based on use):

  - **General use of ripserr:** Wadhwa RR, Piekenbrock M, Scott JG.
    ripserr: Calculate Persistent Homology with Ripser-based Engines;
    version 0.1.0. URL <https://github.com/rrrlw/ripserr>.
  - **Calculation using Vietoris-Rips complex:** Bauer U. Ripser:
    Efficient computation of Vietoris-Rips persistence barcodes. 2019;
    arXiv: 1908.02518.
  - **Calculation using cubical complex:** Kaji S, Sudo T, Ahara K.
    Cubical Ripser: Software for computing persistent homology of image
    and volume data. 2020; arXiv: 2005.12692.

## Contribute

To contribute to ripserr, you can create issues for any bugs/suggestions
on the [issues page](https://github.com/rrrlw/ripserr/issues). You can
also fork the ripserr repository and create pull requests to add useful
features.
