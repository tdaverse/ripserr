# s3 generic
calc_phom <- function(x, ...) {
  UseMethod("calc_phom")
}

# time series objects (class ts, also include others like xts?)
calc_phom.ts <- function() {
  
}

# data frames (also include others like tibbles?)
calc_phom.data.frame <- function() {
  
}

# arrays (also include matrix and others?) - probs not matrix b/c r-core transition
calc_phom.array <- function() {
  
}