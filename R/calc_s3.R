# s3 generic
#' @rdname calc_phom
#' @export calc_phom
calc_phom <- function(x, ...) {
  UseMethod("calc_phom")
}

# assume object of class numeric (single vector) is a time series
# quasi-attractor transform from:
#     Umeda Y. Time Series Classification via Topological Data Analysis.
#     Transactions of the Japanese Society for Artificial Intelligence. 2017;
#     32(3): DG72 1-12. doi: 10.1527/tjsai.D-G72
#' @rdname calc_phom
#' @export calc_phom.numeric
#' @export
calc_phom.numeric <- function(x,
                              data_dim = 2L, max_dim = 1L,
                              initial_term = 1L,
                              smp_lag = 1L) {
  
  seq_index <- seq(from = initial_term,
                   by = smp_lag,
                   to = length(x))
  seq_index <- x[seq_index]
  
  tmp_mat <- matrix(NA,
                    ncol = data_dim,
                    nrow = length(seq_index) - data_dim + 1)
  
  for (curr_index in seq_len(ncol(tmp_mat))) {
    tmp_mat[, curr_index] <- seq_index[curr_index : (curr_index + nrow(tmp_mat) - 1)]
  }
  
  tmp_df <- as.data.frame(tmp_mat)
  
  return(ripserr::calc_phom.data.frame(tmp_df))
}

# time series objects (class ts; include xts, zoo, etc. in TDAstats)
# switch quasi-attractor code from numeric to ts; then add param method "quasi" to ts
#' @rdname calc_phom
#' @export calc_phom.ts
#' @export
calc_phom.ts <- function(x,
                         data_dim = 2L, max_dim = 1L,
                         initial_term = 1L,
                         smp_lag = 1L) {
  # convert ts object to appropriately formatted data frame
  intermediary <- as.numeric(x)
  
  #####INSERT HERE#####
  
  # calculate using calc_phom.data.frame
  return(ripserr::calc_phom.numeric(intermediary,
                                    data_dim = data_dim,
                                    max_dim = max_dim,
                                    initial_term = initial_term,
                                    smp_lag = smp_lag))
}

# for dist objects
# can only be vietoris_rips
#' @rdname calc_phom
#' @export calc_phom.dist
#' @export
calc_phom.dist <- function(x,
                           max_dim = 1L, threshold = -1, p = 2L,
                           standardize = FALSE, return_format = "df") {
  ripserr::vietoris_rips(dataset = x,
                         dim = max_dim,
                         threshold = threshold,
                         p = p,
                         standardize = standardize,
                         return_format = return_format)
}

# data frames (include tibble in TDAstats)
# can only be vietoris_rips (2d cubical should be in array/matrix format)
#' @rdname calc_phom
#' @export calc_phom.data.frame
#' @export
calc_phom.data.frame <- function(x,
                                 max_dim = 1L, threshold = -1, p = 2L,
                                 standardize = FALSE, return_format = "df") {
  ripserr::vietoris_rips(dataset = x,
                         dim = max_dim,
                         threshold = threshold,
                         p = p,
                         standardize = standardize,
                         return_format = return_format)
}

# arrays (includes matrix)
#' @rdname calc_phom
#' @export calc_phom.array
#' @export
calc_phom.array <- function(x,
                            threshold = 9999, method = 0,
                            standardize = FALSE, return_format = "df") {
  ripserr::cubical(dataset = x,
                   threshold = threshold,
                   method = method,
                   standardize = standardize,
                   return_format = return_format)
}