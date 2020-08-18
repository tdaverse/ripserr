# ripserr: Calculate Persistent Homology of Vietoris-Rips and Cubical Complexes using Ripser in R

[![Travis-CI Build Status](https://travis-ci.org/rrrlw/ripserr.svg?branch=master)](https://travis-ci.org/rrrlw/ripserr)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/rrrlw/ripserr?branch=master&svg=true)](https://ci.appveyor.com/project/rrrlw/ripserr)

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN version](http://www.r-pkg.org/badges/version/ripserr)](https://CRAN.R-project.org/package=ripserr)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/ripserr)](https://CRAN.R-project.org/package=ripserr)

## Overview

ripserr ports the [Ripser](https://arxiv.org/1908.02518) and [Cubical Ripser](https://arxiv.org/abs/2005.12692) persistent homology calculation engines from C++ via [Rcpp](https://CRAN.R-project.org/package=Rcpp).
It can be used as a convenient and rapid calculation tool in topological data analysis pipelines.

## Installation

```r
# install development version
devtools::install_github("rrrlw/ripserr")

# install from CRAN (not available yet)
# install.packages("ripserr")
```

## Sample code

## Functionality

1. Calculation of persistent homology of Vietoris-Rips complexes using Ripser.
1. Calculation of persistent homology of cubical complexes using Cubical Ripser.

## Contribute
