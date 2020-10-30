#####AUTO-ERROR MESSAGES#####
error_class <- function(var, param_name, desired_class) {
  if (all(!(desired_class %in% class(var)))) {
    stop(paste0(param_name, " parameter must be of class ", desired_class,
                ", passed value (class) = ", var, " (", class(var), ")"))
  }
}

error_integer <- function(x, param_name) {
  error_class(x, param_name, c("integer", "numeric"))
  
  if (!close_to_integer(x)) {
    stop(paste(param_name, "parameter must be an integer, passed value =",
               x, "with class", class(x)))
  }
}

#####PARAMETER VALIDATION FUNCTIONS#####
# confirm that x is within epsilon distance from an integer
close_to_integer <- function(x, epsilon = 1e-6) {
  return(abs(x - round(x)) < epsilon)
}

# make sure parameters for vietoris_rips make sense
validate_params_vr <- function(max_dim, threshold, p) {
  # stuff for max_dim
  error_integer(max_dim, "max_dim")
  
  if (max_dim < 0) {
    stop(paste("max_dim parameter must be nonnegative, passed value =",
               max_dim))
  }
  
  # stuff for threshold
  error_class(threshold, "threshold", c("integer", "numeric"))
  
  # stuff for p
  # primality is checked in C++
  error_integer(p, "p")
}

# make sure parameters for vietoris_rips time series make sense
validate_params_ts_vr <- function(vec_len,
                                  data_dim, max_dim,
                                  dim_lag, sample_lag,
                                  method) {
  # data_dim checks
  error_integer(data_dim, "data_dim")
  if (data_dim < 2) {
    stop(paste("data_dim parameter must be at least 2, passed value =",
               data_dim))
  }
  
  # max_dim checks
  # error_integer(max_dim, "max_dim")
  # if (max_dim < 0) {
  #   stop(paste("max_dim parameter must be nonnegative, passed value =",
  #              max_dim))
  # } else if (max_dim > data_dim - 1) {
  #   stop(paste("max_dim must be at least 1 less than data_dim, passed",
  #              "data_dim =", data_dim, "and max_dim =", max_dim))
  # }
  
  # dim_lag checks
  error_integer(dim_lag, "dim_lag")
  if (dim_lag < 1) {
    stop(paste("dim_lag must be a positive integer, passed value =", dim_lag))
  } else if (dim_lag > vec_len) {
    stop(paste("dim_lag must be less than the number of elements in the time",
               "series, passed dim_lag =", dim_lag, "and time series length =",
               vec_len))
  }
  
  # sample_lag checks
  error_integer(sample_lag, "sample_lag")
  if (sample_lag < 1) {
    stop(paste("sample_lag must be a positive integer, passed value =",
               sample_lag))
  } else if (sample_lag > vec_len) {
    stop(paste("sample_lag must be less than the number of elements in the",
               "time series, passed dim_lag =", dim_lag, "and time series",
               "length =", vec_len))
  }
  
  # math checks for combos
  n_rows <- ceiling(vec_len / sample_lag) - (data_dim - 1) * dim_lag
  if (n_rows < 3) {
    stop(paste0("With the provided parameters, the number of rows in the ",
         "resulting matrix will be ", n_rows, " (less than 2), passed values ",
         "of relevant variables are: time series length = ", vec_len,
         "; sample_lag = ", sample_lag, "; data_dim = ", data_dim,
         "; dim_lag = ", dim_lag))
  }
  
  # method checks; only quasi-attractor ("qa") implemented so far
  if (!(method %in% c("qa"))) {
    stop(paste(method, "is not a valid option for the method parameter"))
  }
}

validate_mat_vr <- function(dataset) {
  # correct class
  error_class(dataset, "dataset", "matrix")
  
  # min size
  if (nrow(dataset) < 2) {
    stop(paste("dataset parameter must have at least 2 rows, number of rows in",
               "passed dataset =", nrow(dataset)))
  } else if (ncol(dataset) < 2) {
    stop(paste("dataset parameter must have at least 2 columns, number of",
               "columns in passed dataset =", ncol(dataset)))
  }
  
  # no missing elements
  if (!all(stats::complete.cases(dataset))) {
    stop(paste("dataset parameter must not have any missing values, rows with",
               "missing values in passed data frame =",
               which(!stats::complete.cases(dataset))))
  }
  
  # all numeric elements
  if (!is.numeric(dataset)) {
    stop("dataset must contain numeric values, passed dataset has class =",
         class(dataset))
  }
}

validate_dist_vr <- function(dataset) {
  # correct class
  error_class(dataset, "dataset", "dist")
  
  # min size
  if (attr(dataset, "Size") < 3) {
    stop("dataset parameter must contain at least 3 elements, passed dist",
         "object contains", attr(dataset, "Size"), "elements")
  }
  
  # no missing elements
  if (!all(stats::complete.cases(dataset))) {
    stop(paste("dataset parameter must not have any missing values, elements",
               "that are missing values inmissing values in passed dist =",
               which(!stats::complete.cases(dataset))))
  }
  
  # all numeric elements
  if (!all(is.numeric(dataset))) {
    stop(paste("dataset parameter must only contain numeric values, elements",
               "in passed dist object that are not numeric =",
               which(!is.numeric(dataset))))
  }
}

#####DATA FORMATTING#####
# convert Ripser C++/Rcpp output into an appropriate data frame
ripser_vec_to_df <- function(ans_vec) {
  # first convert to matrix (easily handled by `matrix` function)
  ans_mat <- matrix(ans_vec,
                    byrow = TRUE,
                    ncol = 3)
  colnames(ans_mat) <- c("dimension", "birth", "death")
  
  # convert to df format and fix col class
  ans_df <- as.data.frame(ans_mat)
  ans_df$dimension <- as.integer(ans_df$dimension)
  
  # add "PHom" class
  class(ans_df) <- c("PHom", "data.frame")
  
  # finally return properly formatted object
  return(ans_df)
}

# convert numeric vector (time series) to matrix for persistent homology
#   calculation based on quasi-attractor method in:
#     Umeda Y. Time Series Classification via Topological Data Analysis.
#     Transactions of the Japanese Society for Artificial Intelligence. 2017;
#     32(3): DG72 1-12. doi: 10.1527/tjsai.D-G72
numeric_to_quasi_attractor <- function(vec, data_dim,
                                       dim_lag, sample_lag) {
  # sequence of terms to be included in final series
  sample_seq <- seq(from = 1, to = length(vec), by = sample_lag)
  
  # setup matrix
  ans_mat <- matrix(NA,
                    ncol = data_dim,
                    nrow = length(sample_seq) - dim_lag * (data_dim - 1))
  
  # fill in matrix column-wise
  for (curr_col in seq_len(ncol(ans_mat))) {
    start_seq <- 1 + dim_lag * (curr_col - 1)
    end_seq <- start_seq + nrow(ans_mat) - 1
    
    curr_seq <- sample_seq[start_seq : end_seq]
    curr_val <- vec[curr_seq]
    
    ans_mat[, curr_col] <- curr_val
  }
  
  # return
  return(ans_mat)
}
