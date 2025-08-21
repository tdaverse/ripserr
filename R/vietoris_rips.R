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
#' matrix using the sliding window embedding introduced in Perea & Harer (2015)
#' <doi:10.1007/s10208-014-9206-z> and used to obtain quasi-attractors in Umeda
#' (2017) <doi:10.1527/tjsai.D-G72>. Persistent homology of the resulting matrix
#' is then calculated. (NB: If a multi-time series is unclassed, then method
#' dispatch will pass it to `vietoris_rips.matrix`).
#' 

#' @param dataset object on which to calculate persistent homology
#' @param ... other relevant parameters
#' @rdname vietoris_rips
#' @export vietoris_rips
#' @return `PHom` object
#' @examples
#'
#' # create a 2-d point cloud of a circle (100 points)
#' num.pts <- 100
#' rand.angle <- runif(num.pts, 0, 2*pi)
#' pt.cloud <- cbind(cos(rand.angle), sin(rand.angle))
#'
#' # calculate persistent homology (num.pts by 3 numeric matrix)
#' pers.hom <- vietoris_rips(pt.cloud)
#' 
#' # sliding window persistent homology
#' ( ld.phom <- vietoris_rips(ldeaths, data_dim = 12) )
#' plot(
#'   ld.phom$birth, ld.phom$death,
#'   pch = ld.phom$dimension + 1, col = ld.phom$dimension + 1,
#'   xlim = c(0, 3200), ylim = c(0, 3200), asp = 1
#' )
#' abline(a = 0, b = 1)
#' 
#' # multiple time series; window defaults to period
#' ( sb.phom <- vietoris_rips(Seatbelts[, 2:4]) )
#' plot(
#'   sb.phom$birth, sb.phom$death,
#'   pch = sb.phom$dimension + 1, col = sb.phom$dimension + 1,
#'   xlim = c(0, 1600), ylim = c(0, 1600), asp = 1
#' )
#' abline(a = 0, b = 1)
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

#' @param max_dim maximum dimension of persistent homology features to be
#'   calculated
#' @param dim deprecated; passed to `max_dim` or ignored if `max_dim` is
#'   specified
#' @param threshold maximum simplicial complex diameter to explore
#' @param p prime field in which to calculate persistent homology
#' @rdname vietoris_rips
#' @export vietoris_rips.matrix
#' @export
vietoris_rips.matrix <- function(
    dataset,
    max_dim = 1L,
    threshold = -1,
    p = 2L,
    dim = NULL,
    ...
) {
  
  # shortcut for special case (only 1 row should return empty PHom)
  if (nrow(dataset) == 1L) return(new_PHom())
  
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
  ans <- new_PHom(ripser_ans_to_df(ans))
  
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
  ans <- new_PHom(ripser_ans_to_df(ans))
  
  # return
  return(ans)
}

#' @aliases vietoris_rips.numeric vietoris_rips.ts
#' @importFrom stats tsp
#' @param data_dim desired end data dimension (for `"ts"`, defaults to obs/time
#'   if > 1)
#' @param dim_lag time series lag factor between dimensions
#' @param sample_lag time series lag factor between samples (rows)
#' @param method currently only allows `"qa"` (quasi-attractor method)
#' @rdname vietoris_rips
#' @export vietoris_rips.numeric
#' @export
vietoris_rips.numeric <- function(
    dataset,
    data_dim = 2L,
    dim_lag = 1L, sample_lag = 1L,
    method = "qa",
    ...
) {
  dataset <- as.matrix(dataset)
  
  # ensure valid arguments passed
  validate_params_ts_vr(
    vec_len = nrow(dataset),
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
vietoris_rips.ts <- function(dataset, data_dim = max(tsp(dataset)[3], 2), ...) {
  
  # convert to numeric and pass to vietoris_rips.numeric
  ans <- vietoris_rips.numeric(as.matrix(dataset), data_dim = data_dim, ...)
  
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
