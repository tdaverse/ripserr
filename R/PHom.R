# PHom: S3 class for basic persistence data (PHom objects should also have "data.frame" %in% class)
# Structure of data = set of features, each of which has:
# - dimension <int>
# - birth <dbl>
# - death <dbl>

#####BACKEND BASICS#####
new_PHom <- function(x = data.frame(dimension = integer(0),
                                    birth = double(0),
                                    death = double(0)),
                     dim_col = 1,
                     b_col = 2,
                     d_col = 3) {
  # ensure valid parameters (do this in PHom helper instead?)
  
  # construct df for PHom object
  ans <- data.frame(dimension = as.integer(x[, dim_col]),
                    birth = as.double(x[, b_col]),
                    death = as.double(x[, d_col]))
    
  # structure class appropriately
  structure(ans,
            class = c("PHom", "data.frame"))
}

validate_PHom <- function(x, error = TRUE) {
  # complete checks
  
  # return original object if all okay
  x
}

# export helper
PHom <- function(x, dim_col = 1, birth_col = 2, death_col = 3) {
  # any parameter checks needed?
  
  # return PHom object
  validate_PHom(
    new_PHom(x = x,
             dim_col = dim_col,
             b_col = birth_col,
             d_col = death_col)
  )
}

#####CONVERTER/CHECKER#####
as.PHom <- function(x, dim_col = 1, birth_col = 2, death_col = 3) {
  x <- as.data.frame(x)
  
  return(PHom(x))
}

is.PHom <- function(x) {
  # use validate to implement checks
  return(
    validate_PHom(x = x, error = FALSE)
  )
}

#####S3 GENERICS#####
