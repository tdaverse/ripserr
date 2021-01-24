library(RJSONIO)
library(hash)
library(R.utils)
library(rdryad)
library(tools)

brazilStateCodeLookup <- hash(
  "Acre" = "AC",
  "Alagoas" = "AL",
  "Amapá" = "AP",
  "Amazonas" = "AM",
  "Bahia" = "BA",
  "Ceará" = "CE",
  "Distrito Federal" = "DF",
  "Espírito Santo" = "ES",
  "Goiás" = "GO",
  "Maranhão" = "MA",
  "Mato Grosso" = "MT",
  "Mato Grosso do Sul" = "MS",
  "Minas Gerais" = "MG",
  "Pará" = "PA",
  "Paraíba" = "PB",
  "Paraná" = "PR",
  "Pernambuco" = "PE",
  "Piauí" = "PI",
  "Rio de Janeiro" = "RJ",
  "Rio Grande do Norte" = "RN",
  "Rio Grande do Sul" = "RS",
  "Rondônia" = "RO",
  "Roraima" = "RR",
  "Santa Catarina" = "SC",
  "São Paulo" = "SP",
  "Sergipe" = "SE",
  "Tocantins" = "TO"
)

# Set the default sleep delay as required by the Nominatim API documentation to not make more than one request per second.
nominatimSleepDelay = 1.1

#' Pre-pends a message with a timestamp.
#'
#' @param message The message to print.
#'
#' @examples
#' printMessageWithTimestamp("This is my example message.")
printMessageWithTimestamp <- function(message) {
  print(paste0(Sys.time(), " ", message))
}

#' Reverse geocode a coordinate pair.
#'
#' This converts the coordinates into state and country names.
#'
#' @param longitude Longitude of the coordinate. Sometimes labeled as the "X" coordinate.
#' @param latitude Latitude of the coordinate. Sometimes labeled as the "Y" coordinate.
#' @param email Email address of the person running this script. Required by the Nominatim API Terms of Service.
#'
#' @return This will return a list with names "state_name", "state_code", "country_name", and "country_code"
#'
#' @examples
#' reverseGeocode(-32.21, -52.38, "test@@asdf.com")
reverseGeocode <- function (longitude, latitude, email)
{
  # Pause as required by the Nominatim API.
  Sys.sleep(nominatimSleepDelay)
  
  url <- paste0(
    "https://nominatim.openstreetmap.org/reverse?format=json&zoom=5&addressdetails=1&lon=",
    longitude,
    "&lat=",
    latitude,
    "&email=",
    email
  )
  
  # Request URL with user agent settings, as required by the Nominatim API documentation.
  userAgent = paste0("R library ripserr demonstration script. One time geocoding of data requested by ",
                     email)
  withr::with_options(
    list(HTTPUserAgent = userAgent),
    response <- readLines(url, encoding = "UTF-8", warn = FALSE)
  )
  response <- fromJSON(response)
  
  state <-
    tryCatch(
      response$address["state"],
      error = function(e)
        NULL
    )
  country <-
    tryCatch(
      response$address["country"],
      error = function(e)
        NULL
    )
  country_code <-
    tryCatch(
      response$address["country_code"],
      error = function(e)
        NULL
    )
  
  result <- list()
  result["state_name"] <- state
  result["state_code"] <- brazilStateCodeLookup[[state]]
  result["country_name"] <- country
  result["country_code"] <- country_code
  
  return(result)
}


#' Geocode every row in a file.
#'
#' For demonstration purposes here, only rows that have the field "COUNTRY" set
#' to "Brasil" will be geocoded. All others will be ignored.
#'
#' @param inputFilePath CSV file to read lat/lon data from.
#' @param outputFilePath CSV file path to write geocoded results to. If this file exists, it will be deleted.
#' @param email Email address of the person running this script. Required by the Nominatim API Terms of Service.
#'
#' @examples
#' processFile("path/to/input.csv", "path/to/output.csv", "test@@asdf.com")
processFile <- function(inputFilePath, outputFilePath, email)
{
  if (file.exists(outputFilePath)) {
    # Remove past run results
    printMessageWithTimestamp("Warning: deleting previous run results.")
    file.remove(outputFilePath)
  }
  
  lineCount = countLines(inputFilePath)
  estimatedProcessingMinutes = lineCount * nominatimSleepDelay / 60
  printMessageWithTimestamp(
    paste0(
      "Based on your file size, this process is expected to take no more than ",
      round(estimatedProcessingMinutes, 2),
      " minutes. This is likely a significant over-estimate."
    )
  )
  
  header <- read.csv(inputFilePath, nrows = 1, header = FALSE)
  headerNeedsWriting = TRUE
  
  printMessageWithTimestamp("Reverse geocoding has started.")
  for (rowIndex in 1:lineCount) {
    if (0 == rowIndex %% 10) {
      printMessageWithTimestamp(paste0(
        "Progress: ",
        rowIndex,
        " rows of ",
        lineCount,
        " have been processed."
      ))
    }
    
    outputRow <-
      read.csv(
        inputFilePath,
        skip = rowIndex,
        nrows = 1,
        header = FALSE,
        col.names = header
      )
    
    if (!identical(outputRow$COUNTRY, "Brazil")) {
      next
    }
    
    x  <- outputRow$X
    y  <- outputRow$Y
    
    if (identical(x, logical(0)) || identical(y, logical(0))) {
      next
    }
    geocodeResult <- reverseGeocode(x, y, email)
    
    outputRow["GEOCODED_STATE_NAME"] <- geocodeResult$state_name
    outputRow["GEOCODED_STATE_CODE"] <- geocodeResult$state_code
    outputRow["GEOCODED_COUNTRY_NAME"] <- geocodeResult$country_name
    outputRow["GEOCODED_COUNTRY_CODE"] <- geocodeResult$country_code
    
    suppressWarnings(
      write.table(
        outputRow,
        outputFilePath,
        append = TRUE,
        sep = ",",
        col.names = headerNeedsWriting,
        row.names = FALSE,
        quote = FALSE,
        fileEncoding = "UTF8"
      )
    )
    
    headerNeedsWriting = FALSE
  }
  
  printMessageWithTimestamp("Reverse geocoding has completed.")
}

#' Run the script in interactive mode.
#'
#' Interactive mode will prompt the user for input files and email.
executeInteractiveScript <- function() {
  email = readline(
    "Please provide your email. We need this to comply with terms of service for the nominatim geocoder). "
  )
  
  downloadDryadData = readline(
    "Do you want to download data from Dryad Data? Y = download the data directly from datadryad.org, N = enter your own file path. "
  )
  
  inputPath = "data-raw/geocode-data/sample_zika_data.csv"
  
  if ("Y" == downloadDryadData || "y" == downloadDryadData) {
    dryad_doi = "10.5061/dryad.47v3c"
    downloaded_documents = dryad_download(dryad_doi)
    
    printMessageWithTimestamp("Downloaded the following CSV files from dryad:")
    for (filePath in downloaded_documents[[dryad_doi]]) {
      if ("csv" != file_ext(basename(filePath))){
        next
      }
      
      print(paste0("- ", basename(filePath)))
    }
    
    for (filePath in downloaded_documents[[dryad_doi]]) {
      if ("csv" != file_ext(basename(filePath))){
        next
      }
      
      shouldProcess = readline(paste0(
        "Do you want to process this file? Y/n: ",
        basename(filePath),
        " "
      ))
      
      if ("Y" == shouldProcess || "y" == shouldProcess) {
        inputPath = filePath
        
        break
      }
    }
  } else {
    filepath = readline(
      "Please provide the fully-qualified path to your local csv file. Leave blank to use the smaller sample data set. "
    )
    
    if ("" != filepath) {
      inputPath = filepath
    }
  }
  
  outputPath = paste0(tools::file_path_sans_ext(inputPath), "_annotated.csv")
  
  printMessageWithTimestamp(paste0("Processing data in ", inputPath))
  printMessageWithTimestamp(paste0("Writing results to ", outputPath))
  
  processFile(inputPath, outputPath, email)
  
  printMessageWithTimestamp(paste0("Results were written to ", outputPath))
}

executeInteractiveScript()
