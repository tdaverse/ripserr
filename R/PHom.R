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
                 "paseed value has class =", class(x)))
    } else {
      return(FALSE)
    }
  }
  if (ncol(x) != 3) {
    if (error) {
      stop(paste("PHom objects must be data frames with exactly 3 columns;",
                 "passed object has", ncol(x), "columns."))
    } else {
      return(FALSE)
    }
  }
  if (!("integer" %in% class(x$dimension))) {
    if (error) {
      stop(paste("`dimension` column in PHom objects must have class integer;",
                 "`dimension` column in passed object has class =",
                 class(x$dimension)))
    } else {
      return(FALSE)
    }
  }
  if (!("numeric" %in% class(x$birth) & "numeric" %in% class(x$death))) {
    if (error) {
      stop(paste("`birth` and `death` columns in PHom objects must have",
                 "class integer; respective columns in passed object have",
                 "classes =", class(x$birth), "and", class(x$death)))
    } else {
      return(FALSE)
    }
  }
  
  # make sure all deaths are after corresponding births
  if (!all(x$birth < x$death)) {
    if (error) {
      stop(paste("In PHom objects, all births must be before corresponding",
                 "deaths."))
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
                 "provided data frame; passed value =", val))
    }
  } else if (is.character(val)) {
    if (!(val %in% colnames(df))) {
      stop(paste(val_name, "parameter must be a valid column name for the",
                 "provided data frame; passed value =", val))
    }
  } else {
    stop(paste(val_name, "parameter must either be a valid column name or",
               "index; passed value =", val))
  }
}

# export helper
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
# base::print S3 generic
print.PHom <- function(x) {
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
