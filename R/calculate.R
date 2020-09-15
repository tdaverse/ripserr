#' Calculate Persistent Homology of a Point Cloud
#'
#' Calculates the persistent homology of a point cloud, as represented by
#' a Vietoris-Rips complex. This function is an R wrapper for Ulrich Bauer's
#' Ripser C++ library for calculating persistent homology. For more
#' information on the C++ library, see <https://github.com/Ripser/ripser>.
#'
#' The `mat` parameter should be a numeric matrix with each row corresponding
#' to a single point, and each column corresponding to a single dimension. Thus,
#' if `mat` has 50 rows and 5 columns, it represents a point cloud with 50 points
#' in 5 dimensions. The `dim` parameter should be a positive integer.
#' Alternatively, the `mat` parameter could be a distance matrix (upper
#' triangular half is ignored); note: `format` should be specified as "distmat".
#'
#' @param dataset numeric matrix containing point cloud or distance matrix
#' @param dim maximum dimension of features to calculate
#' @param threshold maximum diameter for computation of Vietoris-Rips complexes
#' @param p number of the prime field Z/pZ to compute the homology over
#' @param format  format of `mat`, either `"cloud"` for point cloud or `"distmat"` for distance matrix
#' @param standardize boolean determining whether point cloud size should be standardized
#' @param return_format defaults to `"df"`, returning a data frame;
#'   if `mat`, returns a numeric matrix
#' @return 3-column matrix or data frame, with each row representing a TDA feature
#' @importFrom stats complete.cases
#' @export
#' @examples
#'
#' # create a 2-d point cloud of a circle (100 points)
#' num.pts <- 100
#' rand.angle <- runif(num.pts, 0, 2*pi)
#' pt.cloud <- cbind(cos(rand.angle), sin(rand.angle))
#'
#' # calculate persistent homology (num.pts by 3 numeric matrix)
#' pers.hom <- vietoris_rips(pt.cloud)
vietoris_rips <- function(dataset, dim = 1, threshold = -1, p = 2L, format = "cloud",
                               standardize = FALSE, return_format = "df") {

  # make sure dim is an integer greater than or equal to zero
  if (!is.numeric(dim)) {
    stop(paste("dim needs to be of class integer or numeric; passed class =",
               class(dim)))
  }
  dim <- as.integer(dim)
  if (dim < 0) {
    stop(paste("dim cannot be negative; passed value =",
               dim))
  }

  # make sure threshold is of type numeric
  if (!is.numeric(threshold)) {
    stop(paste("threshold parameter must be of class integer or numeric; current class =",
               class(threshold)))
  }

  # make sure p is integer; whether p is prime is checked in C++
  if (!is.numeric(p)) {
    stop(paste("p parameter must be of class integer or numeric; current class =",
               class(p)))
  }
  p <- as.integer(p)
  
  # make sure standardize is boolean
  if (!is.logical(standardize)) {
    stop(paste("standardize parameter must be of class logical; passed value =",
               standardize))
  }

  ## If dist object passed in directly, use that, otherwise check various conditions that
  ## mat is a valid matrix object
  if (methods::is(dataset, "dist")) {
    format <- "distmat"
    stopifnot(all(is.numeric(dataset)))
    ans_vec <- ripser_cpp_dist(dataset, dim, threshold, p)
  
  ## not a dist object
  } else {
    # coerce mat into matrix to work with class object such as dist class object
    mat <- as.matrix(dataset)

    # make sure matrix has at least 2 columns and at least 2 rows
    if (nrow(mat) < 2 | ncol(mat) < 2) {
      stop("Point cloud must have at least 2 points and at least 2 dimensions.")
    }

    # make sure matrix contains only numeric (includes integer) values
    if (!all(is.numeric(mat))) {
      stop("Point cloud must contain values of class `numeric` or `integer` only.")
    }

    # make sure there are no NAs in matrix
    if (sum(stats::complete.cases(mat)) < nrow(mat)) {
      stop("Point cloud has missing values.")
    }

    # make sure format is either "cloud" or "distmat"
    if (!(format %in% c("cloud", "distmat"))) {
      stop("format parameter should be either \"cloud\" or \"distmat\"")
    }
    int.format <- ifelse(format == "cloud", 0, 1)

    # standardize if necessary
    # method: independently for each coordinate axis: (x - min) / (max - min)
    if (standardize) {
      for (i in 1:ncol(mat)) {
        min.val <- min(mat[, i])
        max.val <- max(mat[, i])

        # skip if only one unique value in this column (all identical)
        if (min.val == max.val) next

        mat[, i] <- (mat[, i] - min.val) / (max.val - min.val)
      }
    }
    # actually do work
    ans_vec <- ripser_cpp(mat, dim, threshold, p, int.format)
  }

  # format properly and return
  ans_mat <- matrix(ans_vec,
                    byrow = TRUE,
                    ncol = 3)
  colnames(ans_mat) <- c("dimension", "birth", "death")

  # return the format requested by user
  ans <- switch(return_format,
                "mat" = ans_mat,
                "df" = {
                  temp <- as.data.frame(ans_mat)
                  temp$dimension <- as.integer(temp$dimension)
                  temp
                })
}

#' Calculate Persistent Homology using a Cubical Complex
#'
#' Calculates the persistent homology of a 2- to 4-dimensional numeric array
#' using a Cubical complex. This function is an R wrapper for Takeki Sudo
#' and Kazushi Ahara's Cubical Ripser C++ library. For more information on
#' the C++ library, see <https://github.com/CubicalRipser>.
#'
#' @param dataset numeric array containing pixel/voxel data
#' @param threshold maximum diameter for computation of Cubical complex
#' @param method defaults to 0 for link join; alternatively, can be 1 for
#'   compute pairs. See original Cubical Ripser code at GitHub user
#'   CubicalRipser for details.
#' @inheritParams vietoris_rips
#' @return 3-column matrix with each row representing a TDA feature
#' @export
#' @examples 
#' 
#' # 2-dim example
#' dataset <- rnorm(10 ^ 2)
#' dim(dataset) <- rep(10, 2)
#' cubical_hom2 <- cubical(dataset)
#' 
#' # 3-dim example
#' dataset <- rnorm(8 ^ 3)
#' dim(dataset) <- rep(8, 3)
#' cubical_hom3 <- cubical(dataset)
#' 
#' # 4-dim example
#' dataset <- rnorm(5 ^ 4)
#' dim(dataset) <- rep(5, 4)
#' cubical_hom4 <- cubical(dataset)
cubical <- function(dataset, threshold = 9999, method = 0,
                    standardize = FALSE, return_format = "df") {
  
  # temp var to store before formatting at end
  ans <- NULL
  
  # make sure threshold is of type numeric
  if (!is.numeric(threshold)) {
    stop(paste("threshold parameter must be of class integer or numeric; current class = ",
               class(threshold)))
  }
  
  # make sure dataset is of type array (also includes matrix)
  #   allow type matrix (second condition) for older R versions
  if (!("array" %in% class(dataset) | "matrix" %in% class(dataset))) {
    stop(paste("Data must be of class array or matrix; passed dataset has class =",
               class(dataset)))
  }
  
  # make sure dataset contains only numeric (includes integer) values
  if (!all(is.numeric(dataset))) {
    stop("Data must contain values of class `numeric` or `integer` only.")
  }
  
  # make sure dataset is not too small
  if (prod(dim(dataset)) == 0) {
    stop(paste("dataset must have at least 1 value; passed dataset contains zero"))
  }
  
  # make sure valid method is picked and is numeric
  if (!is.numeric(method)) {
    stop(paste("method must be numeric; class of passed method =",
               class(method)))
  }
  if (!(method %in% c(0, 1))) {
    stop(paste("method must be either 0 (link find) or 1 (compute pairs); passed method =",
                method))
  }
  
  # make sure standardize is boolean
  if (!is.logical(standardize)) {
    stop(paste("standardize parameter must be of class logical; passed value =",
               standardize))
  }
  
  # make all pixels/voxels in [0, 1]
  if (standardize) {
    min_val <- min(dataset)
    max_val <- max(dataset)
    
    dataset <- (dataset - min_val) / (max_val - min_val)
  }
  
  # dim = 2
  if (length(dim(dataset)) == 2) {
    
    # check size limit
    if (dim(dataset)[1] > 2000 |
        dim(dataset)[2] > 1000) {
      stop(paste("Max size for dim 2 = 2000 x 1000; passed size =",
                 dim(dataset)[1], "x", dim(dataset)[2]))
    }
    
    ans <- cubical_2dim(dataset, threshold, method)
  
  # dim = 3
  } else if (length(dim(dataset)) == 3) {
    
    # check size limit
    if (sum(dim(dataset) < 512) < 3) {
      stop(paste("Max size for dim 3 = 512 x 512 x 512; passed size =",
                 dim(dataset)[1], "x", dim(dataset)[2], "x", dim(dataset)[3]))
    }
    
    temp_mat <- dataset
    dim(temp_mat) <- prod(dim(dataset))
    ans <- cubical_3dim(temp_mat, threshold, method,
                        dim(dataset)[1],
                        dim(dataset)[2],
                        dim(dataset)[3])
    
  # dim = 4
  } else if (length(dim(dataset)) == 4) {
    # check size limit
    if (sum(dim(dataset) < 64) < 4) {
      stop(paste("Max size for dim 4 = 64 x 64 x 64 x 64; passed size =",
                 dim(dataset)[1], "x", dim(dataset)[2], "x", dim(dataset)[3], "x", dim(dataset)[4]))
    }
    
    temp_mat <- dataset
    dim(temp_mat) <- prod(dim(dataset))
    ans <- cubical_4dim(temp_mat, threshold, method,
                        dim(dataset)[1],
                        dim(dataset)[2],
                        dim(dataset)[3],
                        dim(dataset)[4])
    
  # idk what to do with `dataset` for now
  } else {
    stop(paste("Can only handle data with dimension = 2, 3, or 4; passed data dimension =", length(dim(dataset))))
  }
  
  ans <- as.data.frame(ans)
  colnames(ans) <- c("dimension", "birth", "death")
  ans$dimension <- as.integer(ans$dimension)
  ans <- switch(return_format,
                "df" = ans,
                "mat" = as.matrix(ans))
  
  return(ans)
}
