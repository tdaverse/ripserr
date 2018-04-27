# ripserr: An R package for Ripser to calculate Vietoris-Rips persistent homology

## Purpose

A popular method of topological data analysis (TDA) is to analyze the topology of a point cloud by building a Vietoris-Rips complex (sometimes shortened simply to Rips complex) and plotting a persistence diagram or topological barcode for further  analysis. However, the process of building a Vietoris-Rips complex is computationally expensive. Ulrich Bauer's Ripser C++ library (https://github.com/Ripser/ripser) dramatically outpaces its competition (Dionysus, DIPHA, GUDHI, Perseus, PHAT) in both speed and memory, but has not been integrated completely with R. Some attempts (e.g. code by [clarkfitzg](https://gist.github.com/clarkfitzg/5f81592cf851e6cfe14187e4cdede56a) and the [RipserOnR](https://github.com/holt0102/RipserOnR) project) to integrate Ripser for R use do exist, but none (to our knowledge) completely integrate Ripser as a functional R package without any manual setup steps.

The purpose of the `ripserr` R package is to eliminate any manual steps for the integration of Ripser and the R programming language (including installation and parsing output), and work with Ripser as though it were an ordinary CRAN-compatible package.

## Use

The idea behind the `ripserr` package is to mimic the Ripser C++ library without adding unnecessary functionality. As such, only a single user-accessible function (`ripser`) is currently implemented, which reads in a point cloud from a local file and returns a 3-column matrix containing the birth, death, and dimension of each feature from the Vietoris-Rips complex of the provided point cloud. We hope to implement a function that takes in a point cloud variable (rather than one stored in a file) in the near future.

Given that the `ripserr` package only calculates Vietoris-Rips persistent homology, we encourage R users to use the `TDA` package in tandem with it for the visualization and comparison of persistence diagrams and topological barcodes. `TDA` provides an R interface for the GUDHI, Dionysus, and PHAT libraries. The `ripserr` package is meant to be a more efficient replacement to its `ripsDiag` function by incorporating the code from the C++ Ripser library. However, the following functions from `TDA` will be quite useful to `ripserr` users:
* `wasserstein(Diag1, Diag2, p, dimension)`: the `Diag1` and `Diag2` parameters correspond to output `n` by 3 matrices from `ripserr`'s `ripser` function; `p` is a parameter related to the Wasserstein metric; `dimension` specifies the dimension for feature comparison.
* `plot.diagram(x, barcode)`: the `x` parameter corresponds to an output `n` by 3 matrix from `ripserr`'s `ripser` function; when the `barcode` parameter is `FALSE`, `plot.diagram` will plot a persistence diagram; when the `barcode` parameter is `TRUE`, `plot.diagram` will plot a topological barcode.
An example of the complementary use of `ripserr` and `TDA` follows.
```r
# create two 2-d point clouds - one a circle, one with randomly distributed points
num.pts <- 100
rand.angle <- runif(num.pts, 0, 2*pi)
pt.cloud1 <- cbind(cos(rand.angle), sin(rand.angle))
pt.cloud2 <- cbind(rnorm(num.pts), rnorm(num.pts))

# calculate persistent homology of the two point clouds with ripserr
homology1 <- ripserr::ripser(pt.cloud1)
homology2 <- ripserr::ripser(pt.cloud2)

# plot the persistence diagram for the first, topological barcode for the second
TDA::plot.diagram(homology1, barcode = FALSE)
TDA::plot.diagram(homology2, barcode = TRUE)

# compare the two homologies using the Wasserstein metric (Betti-1 features)
# (in this case, since p=1, equivalent to the Earth Mover's Distance)
comparison <- TDA::wasserstein(homology1, homology2, p = 1, dimension = 1)
print(comparison)
```
