#' Calculate Persistent Homology via a Vietoris-Rips Complex
#' 
#' @rdname vietoris_rips
#' @export vietoris_rips
# Notes:
# - figure out format from `dataset`
# - return_format will be "df" (opinionated) w/ additional "PHom" S3 class
# - standardize will be a different method (can be connected w/ magrittr pipe)
# - apart from dataset, only `max_dim = 1L, threshold = -1, p = 2L` are needed
vietoris_rips <- function(dataset, ...) {
  UseMethod("vietoris_rips")
}

#' Vietoris-Rips Persistent Homology of a Data Frame
#' 
#' @importFrom magrittr %>%
#' @rdname vietoris_rips
#' @export vietoris_rips.data.frame
#' @export
# must be in cloud format, distmat should go to vietoris_rips.dist
vietoris_rips.data.frame <- function(dataset, ...) {
  # convert to matrix and pass to vietoris_rips.matrix
  ans <- dataset %>%
    as.matrix() %>%
    vietoris_rips.matrix(...)
  
  # return
  return(ans)
}

#' Vietoris-Rips Persistent Homology of a Matrix
#' 
#' @importFrom magrittr %>%
#' @rdname vietoris_rips
#' @export vietoris_rips.matrix
#' @export
vietoris_rips.matrix <- function(dataset,
                                 max_dim = 1L, threshold = -1, p = 2L,
                                 ...) {
  # ensure valid arguments passed
  validate_params_vr(max_dim = max_dim,
                     threshold = threshold,
                     p = p)
  validate_mat_vr(dataset = dataset)
  
  # calculate persistent homology
  ans <- dataset %>%
    ripser_cpp(max_dim, threshold, p, 0) %>%
    ripser_vec_to_df()
  
  # return
  return(ans)
}

#' Vietoris-Rips Persistent Homology of a \code{dist} Object
#' 
#' @importFrom magrittr %>%
#' @rdname vietoris_rips
#' @export vietoris_rips.dist
#' @export
vietoris_rips.dist <- function(dataset,
                               max_dim = 1L, threshold = -1, p = 2L,
                               ...) {
  # ensure valid arguments passed
  validate_params_vr(max_dim = max_dim,
                     threshold = threshold,
                     p = p)
  validate_dist_vr(dataset = dataset)
  
  # calculate persistent homology
  ans <- dataset %>%
    ripser_cpp_dist(max_dim, threshold, p) %>%
    ripser_vec_to_df()
  
  # return
  return(ans)
}

#' Vietoris-Rips Persistent Homology of Time Series
#' 
#' @aliases vietoris_rips.numeric vietoris_rips.ts
#' @importFrom magrittr %>%
#' @rdname vietoris_rips
#' @export vietoris_rips.numeric
#' @export
vietoris_rips.numeric <- function(dataset,
                                  data_dim = 2L,
                                  dim_lag = 1L, sample_lag = 1L,
                                  method = "qa",
                                  ...) {
  # ensure valid arguments passed
  validate_params_ts_vr(vec_len = length(dataset),
                        data_dim = data_dim,
                        dim_lag = dim_lag,
                        sample_lag = sample_lag,
                        method = method)
  
  # construct appropriate matrix from numeric vector (time series)
  converted <- switch(method,
                      qa = numeric_to_quasi_attractor(dataset, data_dim,
                                         dim_lag, sample_lag),
                      stop(paste("invalid method; this line of code should",
                                 "never be reached")))
  
  # pass matrix to vietoris_rips.matrix
  ans <- vietoris_rips.matrix(converted, ...)
  
  # return
  return(ans)
}

#' @importFrom magrittr %>%
#' @rdname vietoris_rips
#' @export vietoris_rips.ts
#' @export
vietoris_rips.ts <- function(dataset, ...) {
  # convert to numeric and pass to vietoris_rips.numeric
  ans <- dataset %>%
    as.numeric() %>%
    vietoris_rips.numeric(...)
  
  # return
  return(ans)
}

#' Default Vietoris-Rips Persistent Homology
#' 
#' @rdname vietoris_rips
#' @export vietoris_rips.default
#' @export
vietoris_rips.default <- function(dataset, ...) {
  stop(paste("no method for the generic vietoris_rips function exists for",
             "objects with class:", class(dataset)))
}