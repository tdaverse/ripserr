#' Calculate Persistent Homology of a Point Cloud
#'
#' Calculates the persistent homology of a point cloud, as represented by
#' a Vietoris-Rips complex. This function is an R wrapper for Ulrich Bauer's
#' Ripser C++ library for calculating persistent homology. For more
#' information on the C++ library, see <https://github.com/Ripser/ripser>.
#'
#' The `mat` parameter should be a numeric matrix with each row corresponding
#' to a single point, and each column corresponding to a single dimension. Thus,
#' if `mat` has 50 rows and 3 columns, it represents a point cloud with 50 points
#' in 3 dimensions. The `dim` parameter should be a positive integer.
#' Alternatively, the `mat` parameter could be a distance matrix (upper
#' triangular half is ignored); note: `format` should be specified as "ldm".
#'
#' @param mat numeric matrix containing point cloud or distance matrix
#' @param dim maximum dimension of features to calculate
#' @param threshold maximum diameter for computation of Vietoris-Rips complexes
#' @param format  format of `mat`, either "cloud" for point cloud or "distmat" for distance matrix
#' @param retval defaults to "mat" (matrix); user can also choose "df" to
#'   instead return a data frame
#' @param use_coef logical determining whether support for coefficients in a
#'   prime field will be enabled
#' @return 3-column matrix or data frame, with each row representing a feature
#' @export
#' @examples
#'
#' # create a 2-d point cloud of a circle (100 points)
#' num_pts <- 100
#' rand_angle <- runif(num_pts, 0, 2 * pi)
#' pt_cloud <- cbind(cos(rand_angle), sin(rand_angle))
#'
#' # calculate persistent homology (num_pts by 3 numeric matrix)
#' ripser(pt_cloud)
ripser <- function(mat,
                   dim = 1,
                   threshold = -1,
                   format = "cloud",
                   retval = "mat",
                   use_coef = FALSE) {

  # make sure matrix has at least 2 columns and at least 2 rows
  if (nrow(mat) < 2 | ncol(mat) < 2) {
    stop("Point cloud must have at least 2 points and at least 2 dimensions.")
  }

  # make sure matrix contains numeric (or integer) values
  # assumption: matrix can only hold objects of one class so only need to check
  #   one element
  temp <- mat[1, 1]
  if (class(temp) != "numeric" && class(temp) != "integer") {
    stop("Point cloud must contain values of class `numeric` or `integer` only.")
  }

  # make sure there are no NAs in matrix
  if (sum(is.na(mat)) > 0) {
    stop("Point cloud has missing values.")
  }

  # actually do work
  ans_vec <- ripser_cpp(mat)

  # format properly and return
  ans_mat <- matrix(ans_vec,
                    byrow = TRUE,
                    ncol = 3)
  colnames(ans_mat) <- c("dimension", "birth", "death")
  return(ans_mat)
}
