# Notes:
# - figure out format from `dataset`
# - return_format will be "df" (opinionated) w/ additional "PHom" S3 class
# - standardize will be a different method (can be connected w/ magrittr pipe)
# - apart from dataset, only `max_dim = 1L, threshold = -1, p = 2L` are needed
cubical <- function(dataset, ...) {
  UseMethod("cubical")
}