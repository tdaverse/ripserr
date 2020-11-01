#' This function is an R wrapper for the CubicalRipser C++ library to calculate
#' persistent homology. For more information on the C++ library, see
#' <https://github.com/CubicalRipser>. For more information on how objects of
#' different classes are evaluated by `cubical`, read the Details section
#' below.
#' 
#' `cubical.array` assumes `dataset` is a lattice, with each element containing
#' the value of the lattice at the point represented by the indices of the
#' element in the `array`.
#' 
#' `cubical.matrix` is redundant for versions of `R` at or after 4.0. For
#' previous versions of `R`, in which objects with class `matrix` do not
#' necessarily also have class `array`, `dataset` is converted to an `array`
#' and persistent homology is then calculated using `cubical.array`.
#'
#' @title Calculating Persistent Homology via a Cubical Complex
#' @param dataset object on which to calculate persistent homology
#' @param ... other relevant parameters
#' @rdname cubical
#' @export cubical
#' @return data frame (also of class `PHom`) with 3 columns and `n` rows, where
#'   column 1 contains feature dimension, column 2 contains feature birth, and
#'   column 3 contains feature death; each row contains 1 of `n` features
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
# Notes:
# - figure out format from `dataset`
# - return_format will be "df" (opinionated) w/ additional "PHom" S3 class
# - standardize will be a different method (can be connected w/ magrittr pipe)
# - apart from dataset, only `threshold = 9999, method = "lj"` are needed
cubical <- function(dataset, ...) {
  UseMethod("cubical")
}

#' @rdname cubical
#' @param threshold maximum simplicial complex diameter to explore
#' @param method either `"lj"` (for Link Join) or `"cp"` (for Compute Pairs);
#'   see Kaji et al. (2020) <arXiv:2005.12692> for details
#' @export cubical.array
#' @export
cubical.array <- function(dataset, threshold = 9999, method = "lj", ...) {
  # ensure valid arguments passed
  validate_params_cub(threshold = threshold,
                      method = method)
  validate_arr_cub(dataset)
  
  # transform method parameter for C++ function
  method_int <- switch(method,
                       lj = 0,
                       cp = 1)
  
  # calculate persistent homology based on dimension of dataset
  ans <- switch(length(dim(dataset)) - 1, # goes from {2,3,4} to {1,2,3} for switch
                # 2-dimensional array
                {
                  cubical_2dim(dataset, threshold, method_int)
                },
                # 3-dimensional array
                {
                  temp_mat <- dataset
                  dim(temp_mat) <- prod(dim(dataset))
                  
                  cubical_3dim(temp_mat, threshold, method_int,
                               dim(dataset)[1],
                               dim(dataset)[2],
                               dim(dataset)[3])
                },
                # 4-dimensional array
                {
                  temp_mat <- dataset
                  dim(temp_mat) <- prod(dim(dataset))
                  
                  cubical_4dim(temp_mat, threshold, method_int,
                               dim(dataset)[1],
                               dim(dataset)[2],
                               dim(dataset)[3],
                               dim(dataset)[4])
                })
  
  # properly format persistent homology output
  ans <- as.data.frame(ans)
  colnames(ans) <- c("dimension", "birth", "death")
  ans$dimension <- as.integer(ans$dimension)
  
  # INEFFICIENT COPYING STEP, TRY TO FIND A WAY AROUND THIS, IF POSSIBLE
  # remove unnecessary feature (dim = -1, birth = min value, death = threshold)
  remove_row <- which(ans$dimension == -1 &
                      close_numeric(ans$birth, min(dataset)) &
                      close_numeric(ans$death, threshold))
  if (is.integer(remove_row) & length(remove_row) == 1) {
    ans <- ans[-remove_row, ]
  }
  
  # add PHom class to data frame
  class(ans) <- c("PHom", "data.frame")
  
  # return
  return(ans)
}

#' @rdname cubical
#' @export cubical.matrix
#' @export
cubical.matrix <- function(dataset, ...) {
  # older R versions (R < 4), where matrix should be converted to array
  #   don't have to do anything for R >= 4
  if (!("array" %in% class(dataset))) {
    dataset <- as.array(dataset,
                        dim = dim(dataset))
  }
  
  # calculate persistent homology using cubical.array
  ans <- cubical.array(dataset, ...)
  
  # return
  return(ans)
}