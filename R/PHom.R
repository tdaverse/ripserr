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
  # parameter validity checked in PHom helper
  
  # construct df for PHom object
  ans <- data.frame(dimension = as.integer(x[, dim_col]),
                    birth = as.double(x[, b_col]),
                    death = as.double(x[, d_col]))
    
  # structure class appropriately
  structure(ans,
            class = c("PHom", "data.frame"))
}

validate_PHom <- function(x, error = TRUE) {
  # check structure of x
  if (!("PHom" %in% class(x) & "data.frame" %in% class(x))) {
    if (error) {
      stop(paste("PHom objects must have classes `PHom` and `data.frame`;",
                 "paseed value has class =", class(x)),
           call. = FALSE)
    } else {
      return(FALSE)
    }
  }
  if (ncol(x) != 3) {
    if (error) {
      stop(paste("PHom objects must be data frames with exactly 3 columns;",
                 "passed object has", ncol(x), "columns."),
           call. = FALSE)
    } else {
      return(FALSE)
    }
  }
  if (!("integer" %in% class(x$dimension))) {
    if (error) {
      stop(paste("`dimension` column in PHom objects must have class integer;",
                 "`dimension` column in passed object has class =",
                 class(x$dimension)),
           call. = FALSE)
    } else {
      return(FALSE)
    }
  }
  if (!("numeric" %in% class(x$birth) & "numeric" %in% class(x$death))) {
    if (error) {
      stop(paste("`birth` and `death` columns in PHom objects must have",
                 "class integer; respective columns in passed object have",
                 "classes =", class(x$birth), "and", class(x$death)),
           call. = FALSE)
    } else {
      return(FALSE)
    }
  }
  
  # make sure all deaths are after corresponding births
  if (!all(x$birth < x$death)) {
    if (error) {
      stop(paste("In PHom objects, all births must be before corresponding",
                 "deaths."),
           call. = FALSE)
    } else {
      return(FALSE)
    }
  }
  
  # return original object if all okay (or TRUE)
  if (error) {
    x
  } else {
    return(TRUE)
  }
}

# validation helper: check if valid column in df is provided
valid_colval <- function(df, val, val_name) {
  if (is.numeric(val)) {
    val <- as.integer(val)
    
    if (!(val %in% seq_len(ncol(df)))) {
      stop(paste(val_name, "parameter must be a valid column index for the",
                 "provided data frame; passed value =", val),
           call. = FALSE)
    }
  } else if (is.character(val)) {
    if (!(val %in% colnames(df))) {
      stop(paste(val_name, "parameter must be a valid column name for the",
                 "provided data frame; passed value =", val),
           call. = FALSE)
    }
  } else {
    stop(paste(val_name, "parameter must either be a valid column name or",
               "index; passed value =", val),
         call. = FALSE)
  }
}

# export helper
#' Persistence Data Container
#' 
#' PHom() creates instances of `PHom` objects, which are convenient containers
#' for persistence data. Generally, data frame (or similar) objects are used
#' to create `PHom` instances with users specifying which columns contain
#' dimension, birth, and death details for each feature.
#' 
#' @param x object used to create `PHom` instance
#' @param dim_col either `integer` representing column index for feature
#'   dimension data or `character` representing column name
#' @param birth_col either `integer` representing column index for feature
#'   birth data or `character` representing column name
#' @param death_col either `integer` representing column index for feature
#'   death data or `character` representing column name
#' @return `PHom` instance
#' @export
#' @examples 
#' # construct data frame with valid persistence data
#' df <- data.frame(dimension = c(0, 0, 1, 1, 1, 2),
#'                  birth = rnorm(6),
#'                  death = rnorm(6, mean = 15))
#' 
#' # create `PHom` instance and print
#' df_phom <- PHom(df)
#' df_phom
#' 
#' # print feature details to confirm accuracy
#' print.data.frame(df_phom)
PHom <- function(x, dim_col = 1, birth_col = 2, death_col = 3) {
  ## basic parameter checks (column nums/names are valid, etc.)
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  
  # check dim_col
  valid_colval(df = x, val = dim_col, val_name = "dim_col")
  
  # check birth_col
  valid_colval(df = x, val = birth_col, val_name = "birth_col")
  
  # check death_col
  valid_colval(df = x, val = death_col, val_name = "death_col")
  
  ## return PHom object
  validate_PHom(
    new_PHom(x = x,
             dim_col = dim_col,
             b_col = birth_col,
             d_col = death_col)
  )
}

#####CONVERTER/CHECKER#####
#' Convert to PHom Object
#' 
#' Converts valid objects to `PHom` instances.
#' 
#' @param x object being converted to `PHom` instance
#' @inheritParams PHom
#' @return `PHom` instance
#' @export
#' @examples 
#' # construct data frame with valid persistence data
#' df <- data.frame(dimension = c(0, 0, 1, 1, 1, 2),
#'                  birth = rnorm(6),
#'                  death = rnorm(6, mean = 15))
#' 
#' # convert to `PHom` instance and print
#' df_phom <- as.PHom(df)
#' df_phom
#' 
#' # print feature details to confirm accuracy
#' print.data.frame(df_phom)
as.PHom <- function(x, dim_col = 1, birth_col = 2, death_col = 3) {
  x <- as.data.frame(x)
  
  return(PHom(x))
}

#' Check PHom Object
#' 
#' Tests if objects are valid `PHom` instances.
#' 
#' @param x object whose `PHom`-ness is being tested
#' @return `TRUE` if `x` is a valid `PHom` object; `FALSE` otherwise
#' @export
#' @examples
#' # create sample persistence data
#' df <- data.frame(dimension = c(0, 0, 1, 1, 1, 2),
#'                  birth = rnorm(6),
#'                  death = rnorm(6, mean = 15))
#' df <- as.PHom(df)
#' 
#' # confirm that persistence data is valid
#' is.PHom(df)
#' 
#' # mess up df object (feature birth cannot be after death)
#' df$birth[1] <- rnorm(1, mean = 50)
#' 
#' # confirm that persistence data is NOT valid
#' is.PHom(df)
is.PHom <- function(x) {
  # use validate to implement checks
  return(
    validate_PHom(x = x, error = FALSE)
  )
}

#####S3 GENERICS#####
# base::print S3 generic
#' Printing Persistence Data
#' 
#' Print a PHom object.
#' 
#' @param x object of class `PHom`
#' @param ... other parameters; ignored
#' @export
#' @examples
#' # create circle dataset
#' angles <- runif(25, 0, 2 * pi)
#' circle <- cbind(cos(angles), sin(angles))
#' 
#' # calculate persistent homology
#' circle_phom <- vietoris_rips(circle)
#' 
#' # print persistence data
#' print(circle_phom)
print.PHom <- function(x, ...) {
  # make sure `x` is a valid PHom object - too expensive for each print??
  stopifnot(is.PHom(x))
  
  dim_counts <- table(x$dimension)
  
  ans1 <- paste("PHom object containing persistence data for",
               nrow(x), "features.")
  
  ans2 <- "Contains:"
  for (i in sort(unique(x$dimension))) {
    curr <- paste0("* ", dim_counts[as.character(i)], " ", i,
                   "-dim feature")
    if (dim_counts[as.character(i)] != 1) {
      curr <- paste0(curr, "s")
    }
    
    ans2 <- paste(ans2, curr, sep = "\n")
  }
  
  ans3 <- paste0("Radius/diameter: min = ", signif(min(x$birth), digits = 5),
                 "; max = ", signif(max(x$death), digits = 5), ".")
  
  cat(paste(ans1, ans2, ans3, sep = "\n\n"))
}

#' First Part of PHom Object
#' 
#' Returns the first part of a `PHom` instance.
#' 
#' @param x object of class `PHom`
#' @param ... other parameters
#' @export
#' @importFrom utils head
#' @examples
#' # create sample persistence data
#' df <- data.frame(dimension = c(0, 0, 1, 1, 1, 2),
#'                  birth = rnorm(6),
#'                  death = rnorm(6, mean = 15))
#' df_phom <- as.PHom(df)
#' 
#' # look at first 3 features
#' head(df_phom)
#' 
#' # look at last 3 features
#' tail(df_phom)
head.PHom <- function(x, ...) {
  x <- as.data.frame(x)
  head(x, ...)
}

#' Last Part of PHom Object
#' 
#' Returns the last part of a `PHom` instance.
#' 
#' @inheritParams head.PHom
#' @export
#' @importFrom utils tail
#' @examples
#' # create sample persistence data
#' df <- data.frame(dimension = c(0, 0, 1, 1, 1, 2),
#'                  birth = rnorm(6),
#'                  death = rnorm(6, mean = 15))
#' df_phom <- as.PHom(df)
#' 
#' # look at first 3 features
#' head(df_phom)
#' 
#' # look at last 3 features
#' tail(df_phom)
tail.PHom <- function(x, ...) {
  x <- as.data.frame(x)
  tail(x, ...)
}
