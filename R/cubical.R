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
#' Note that the superlevel set filtration (`sublevel = FALSE`) will yield a
#' persistence diagram in which `death` precedes `birth`.
#' 
#' @title Calculating Persistent Homology via a Cubical Complex
#' @param dataset object on which to calculate persistent homology
#' @param ... other relevant parameters
#' @rdname cubical
#' @export cubical
#' @return `PHom` object
#' @examples 
#' 
#' # 1-dim example
#' dataset <- rnorm(1500)
#' dim(dataset) <- 1500
#' cubical_hom1 <- cubical(dataset)
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
#' 
#' # sublevel versus superlevel
#' cubical(volcano)
#' cubical(volcano, sublevel = FALSE)
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
#'   see Kaji et al. (2020) <https://arxiv.org/abs/2005.12692> for details
#' @param sublevel logical; whether to take the sublevel set filtration or else
#'   the superlevel set filtration
#' @export cubical.array
#' @export
cubical.array <- function(
    dataset,
    threshold = 9999, method = "lj",
    sublevel = TRUE,
    ...
) {
  # do this before checks since it modifies `dataset`
  if (! is.logical(sublevel) || is.na(sublevel))
    stop("`sublevel` must be `TRUE` or `FALSE`.")
  if (! sublevel) dataset <- -dataset
  
  # ensure valid arguments passed
  validate_params_cub(threshold = threshold,
                      method = method)
  validate_arr_cub(dataset)
  
  # if dataset is 1-dimensional, treat it as 2-dimensional
  if (length(dim(dataset)) == 1L) dim(dataset) <- c(dim(dataset), 1L)
  
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
  
  # INEFFICIENT COPYING STEP, TRY TO FIND A WAY AROUND THIS, IF POSSIBLE
  # remove unnecessary feature (dim = -1, birth = min value, death = threshold)
  remove_row <- which(ans[, 1L] == -1 &
            close_numeric(ans[, 2L], min(dataset)) &
            close_numeric(ans[, 3L], threshold))
  if (is.integer(remove_row) & length(remove_row) == 1L)
    ans <- ans[-remove_row, , drop = FALSE]
  
  if (! sublevel) ans[, c(2L, 3L)] <- -ans[, c(2L, 3L)]
  
  # properly format persistent homology output
  ans <- as.data.frame(ans)
  colnames(ans) <- c("dimension", "birth", "death")
  ans$dimension <- as.integer(ans$dimension)
  
  # convert data frame to a PHom object
  ans <- new_PHom(ans)
  
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

#' @rdname cubical
#' @export cubical.numeric
#' @export
cubical.numeric <- function(dataset, ...) {
  # convert the numeric vector to a 1-dimensional array
  dataset <- as.array(dataset, dim = 1L)
  
  # calculate persistent homology using cubical.array
  ans <- cubical.array(dataset, ...)
  
  return(ans)
}
