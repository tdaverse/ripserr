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
  # ensure valid parameters
  
  # construct df for PHom object
  ans <- data.frame(dimension = as.integer(x[, dim_col]),
                    birth = as.double(x[, b_col]),
                    death = as.double(x[, d_col]))
    
  # structure class appropriately
  structure(ans,
            class = c("PHom", "data.frame"))
}

validate_PHom <- function(x) {
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

#####S3 GENERICS#####
