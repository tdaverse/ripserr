# s3 generic
#' @export
calc_phom <- function(x, ...) {
  UseMethod("calc_phom")
}

# assume object of class numeric (single vector) is a time series
# quasi-attractor transform from:
#     Umeda Y. Time Series Classification via Topological Data Analysis.
#     Transactions of the Japanese Society for Artificial Intelligence. 2017;
#     32(3): DG72 1-12. doi: 10.1527/tjsai.D-G72
#' @export
calc_phom.numeric <- function(x,
                              data_dim = 2L, max_dim = 1L,
                              initial_term = 1L,
                              smp_lag = 1) {
  
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
  
  return(calc_phom.data.frame(tmp_df))
}

# time series objects (class ts; include xts, zoo, etc. in TDAstats)
calc_phom.ts <- function(x,
                         data_dim = 2L, max_dim = 1L,
                         lags = base::seq_len(data_dim - 1)) {
  # convert ts object to appropriately formatted data frame
  intermediary <- data.frame()
  
  # calculate using calc_phom.data.frame
  return(calc_phom.data.frame(intermediary))
}

# data frames (include tibble in TDAstats)
# can only be vietoris_rips (2d cubical should be in array/matrix format)
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
calc_phom.array <- function() {
  
}