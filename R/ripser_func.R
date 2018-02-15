#' Calculate Persistent Homology of a Point Cloud
#'
#' The `mat` parameter should be a numeric matrix with each row corresponding
#' to a single point, and each column corresponding to a single dimension. Thus,
#' if `mat` has 50 rows and 5 columns, it represents a point cloud with 50 points
#' in 5 dimensions.
#'
#' @param mat numeric matrix containing point cloud
#' @export
ripser <- function(mat) {
  ans_mat <- ripser_cpp(mat)
  return(matrix(ans_mat, byrow = TRUE, ncol = 3))
}
