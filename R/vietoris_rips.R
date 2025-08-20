#' @title Calculate Persistent Homology via a Vietoris-Rips Complex
#'
#' @description This function is an R wrapper for the Ripser C++ library to
#'   calculate persistent homology. For more information on the C++ library, see
#'   <https://github.com/Ripser/ripser>. For more information on how objects of
#'   different classes are evaluated by `vietoris_rips`, read the Details
#'   section below.
#'
#'
#' @details
#'
#' `vietoris_rips.data.frame` assumes `dataset` is a point cloud, with each row
#' representing a point and each column representing a dimension.
#'
#' `vietoris_rips.matrix` currently assumes `dataset` is a point cloud (similar
#' to `vietoris_rips.data.frame`). Currently in the process of adding network
#' representation to this method.
#'
#' `vietoris_rips.dist` takes a `dist` object and calculates persistent homology
#' based on pairwise distances. The `dist` object could have been calculated
#' from a point cloud, network, or any object containing elements from a finite
#' metric space.
#'
#' `vietoris_rips.numeric` and `vietoris_rips.ts` both calculate persistent
#' homology of a time series object. The time series object is converted to a
#' matrix using the quasi-attractor method detailed in Umeda (2017)
#' <doi:10.1527/tjsai.D-G72>. Persistent homology of the resulting matrix is
#' then calculated.
#'
#' @importFrom phutil as_persistence
#' @param dataset object on which to calculate persistent homology
#' @param ... other relevant parameters
#' @rdname vietoris_rips
#' @export vietoris_rips
#' @return `"PHom"` or `"persistence"` object
#' @examples
#'
#' # create a 2-d point cloud of a circle (100 points)
#' num.pts <- 100
#' rand.angle <- runif(num.pts, 0, 2*pi)
#' pt.cloud <- cbind(cos(rand.angle), sin(rand.angle))
#'
#' # calculate persistent homology (num.pts by 3 numeric matrix)
#' ( pers.hom <- vietoris_rips(pt.cloud) )
# Notes:
# - figure out format from `dataset`
# - return_format will be "df" (opinionated) w/ additional "PHom" S3 class
# - standardize will be a different method (can be connected w/ magrittr pipe)
# - apart from dataset, only `max_dim = 1L, threshold = -1, p = 2L` are needed
vietoris_rips <- function(dataset, ...) {
  UseMethod("vietoris_rips")
}

#' @rdname vietoris_rips
#' @export vietoris_rips.data.frame
#' @export
vietoris_rips.data.frame <- function(dataset, ...) {
  
  # convert to matrix and pass to vietoris_rips.matrix
  ans <- vietoris_rips.matrix(as.matrix(dataset), ...)
  
  # return
  return(ans)
}

#' @rdname vietoris_rips
#' @param max_dim maximum dimension of persistent homology features to be
#'   calculated
#' @param dim deprecated; passed to `max_dim` or ignored if `max_dim` is
#'   specified
#' @param threshold maximum simplicial complex diameter to explore
#' @param p prime field in which to calculate persistent homology
#' @param return_class class of output object; either `"PHom"` (default; legacy)
#'   or `"persistence"` (from the
#'   **[phutil](https://cran.r-project.org/package=phutil)** package)
#' @export vietoris_rips.matrix
#' @export
vietoris_rips.matrix <- function(
    dataset,
    max_dim = 1L,
    threshold = -1,
    p = 2L,
    dim = NULL,
    return_class = c("PHom", "persistence"),
    ...
) {
  
  # shortcut for special case (only 1 row should return empty PHom)
  if (nrow(dataset) == 1L) {
    return(switch(
      match.arg(return_class),
      PHom = new_PHom(),
      persistence = as_persistence(matrix(NA_real_, nrow = 0L, ncol = 3L))
    ))
  }
  
  # handle `dim` if passed
  if (! is.null(dim)) {
    max_dim_use <- "max_dim" %in% names(match.call())
    warning("`dim` parameter has been deprecated; ",
            if (max_dim_use) "using" else "use",
            " `max_dim` instead.",
            immediate. = TRUE, call. = TRUE)
    if (! max_dim_use) max_dim <- dim
  }
  
  # ensure valid arguments passed
  validate_params_vr(
    max_dim = max_dim,
    threshold = threshold,
    p = p
  )
  validate_mat_vr(dataset = dataset)
  
  # convert no-threshold value
  if (threshold == -1) threshold <- Inf
  
  # convert distance matrix
  dataset <- stats::dist(dataset)
  
  # calculate persistent homology
  ans <- ripser_cpp_dist(dataset, max_dim, threshold, 1., p)
  
  # coerce to 'PHom' class
  ans <- switch(
    match.arg(return_class),
    PHom = new_PHom(ripser_ans_to_df(ans)),
    persistence = as_persistence(ans)
  )
  
  # return
  return(ans)
}

#' @rdname vietoris_rips
#' @export vietoris_rips.dist
#' @export
vietoris_rips.dist <- function(
    dataset,
    max_dim = 1L,
    threshold = -1,
    p = 2L,
    dim = NULL,
    return_class = c("PHom", "persistence"),
    ...
) {
  
  # handle `dim` if passed
  if (! is.null(dim)) {
    max_dim_use <- "max_dim" %in% names(match.call())
    warning("`dim` parameter has been deprecated; ",
            if (max_dim_use) "using" else "use",
            " `max_dim` instead.",
            immediate. = TRUE, call. = TRUE)
    if (! max_dim_use) max_dim <- dim
  }
  
  # ensure valid arguments passed
  validate_params_vr(
    max_dim = max_dim,
    threshold = threshold,
    p = p
  )
  validate_dist_vr(dataset = dataset)
  
  # convert no-threshold value
  if (threshold == -1) threshold <- Inf
  
  # convert distance matrix
  dataset <- dataset
  
  # calculate persistent homology
  ans <- ripser_cpp_dist(dataset, max_dim, threshold, 1., p)
  
  # coerce to 'PHom' class
  ans <- switch(
    match.arg(return_class),
    PHom = new_PHom(ripser_ans_to_df(ans)),
    persistence = as_persistence(ans)
  )
  
  # return
  return(ans)
}

#' @rdname vietoris_rips
#' @param data_dim desired end data dimension
#' @param dim_lag time series lag factor between dimensions
#' @param sample_lag time series lag factor between samples (rows)
#' @param method currently only allows `"qa"` (quasi-attractor method)
#' @export vietoris_rips.numeric
#' @export
vietoris_rips.numeric <- function(
    dataset,
    data_dim = 2L,
    dim_lag = 1L, sample_lag = 1L,
    method = "qa",
    ...
) {
  
  # ensure valid arguments passed
  validate_params_ts_vr(
    vec_len = length(dataset),
    data_dim = data_dim,
    dim_lag = dim_lag,
    sample_lag = sample_lag,
    method = method
  )
  
  # convert numeric vector (as time series) to appropriate matrix
  converted <- switch(
    method,
    qa = numeric_to_quasi_attractor(
      dataset, data_dim,
      dim_lag, sample_lag
    ),
    stop(paste("invalid method; this line of code should never be reached"))
  )
  
  # pass matrix to vietoris_rips.matrix
  ans <- vietoris_rips.matrix(converted, ...)
  
  # return
  return(ans)
}

#' @rdname vietoris_rips
#' @export vietoris_rips.ts
#' @export
vietoris_rips.ts <- function(dataset, ...) {
  
  # convert to numeric and pass to vietoris_rips.numeric
  ans <- vietoris_rips.numeric(as.numeric(dataset), ...)
  
  # return
  return(ans)
}

#' @rdname vietoris_rips
#' @export vietoris_rips.default
#' @export
vietoris_rips.default <- function(dataset, ...) {
  stop(paste("no method for the generic vietoris_rips function exists for",
             "objects with class:", class(dataset)))
}
