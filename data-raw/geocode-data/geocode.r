library(RJSONIO)
library(hash)
library(R.utils)

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
  defaultdownloadurl = "https://uc3-s3mrt1001-prd.s3.us-west-2.amazonaws.com/3877a3a1-b8e2-4c47-9a55-a03572325842/data?response-content-disposition=attachment%3b%20filename%3ddoi_10.5061_dryad.47v3c__v2.zip&response-content-type=application%2fzip&x-amz-security-token=iqojb3jpz2lux2vjeakacxvzlxdlc3qtmijgmeqciblms%2biaicmw%2bvtf3t1ywxfkrewykec%2b6csg2qd7rfq4aibssye%2fqdipxrayzrbfrfzaaeofoiovyiclxhfldvmoncq9awji%2f%2f%2f%2f%2f%2f%2f%2f%2f%2f8beaaaddq1mtgynjkxnde1nyimm3ggolpzburlhpuckpedmwo5kddo1adxqnytewmfuz2svn0redp927zd30da2ec3zba9avwi33%2fwxr2arvd4vwn%2be%2bibfv9o6hhmpzvl91us9o2j6txctic1ok%2fqax%2bnzu1nc1tgmsxr6rqj5aaerpqizlx5my1fqclfldw%2bbi035meqhcjmnv%2fykklr0itxx%2fsan0jqi55vl2u2lh1b%2fyyplx5rd5hdi1%2brxmgqjckuhelacdt2yshhsteynz%2buqoox6wtdz2qh%2f%2fvk3vsorll%2blbfz9ktncy0vvhalr4yvhdy8pvyvwrz415idy7qp0wppkioq1k%2bj5p9qc8stfo53sivzlzxxdofdy8ir84y0qz9dfel4avimm2oqht8aaubhvfjkzsh0mfefqkp%2bvdfebj1wmrkkio4ng8lczrwbngnksxzewu1fopujnqdhashzpwtftan9vy7v494fchqu9g%2flj94ndvpnkxa8ls1gqhct1tawhzbuilwm9abwjedeekprswtpnbnwdhnf61kj3ysrfkioozuqvrimiqkwjjgzgay67ahcocftdmsye3fdn9fixittzp6v%2bydy4bzju63gncyrlei1vwr%2b5ispdemztaudzge3lxsse%2bclaeic3izgn6qbn16fplgi%2bfwb9lold%2bglcecyh0wice7mqqns5lwkgbnqvztciaphpkkzxumpemhzkduz97mgds5yip%2f2b2yt3xklvuedju9b%2bka24n2q%2bwihtsrrbd8cq738ev6op12wkygmbe7r7wmhywpgcm1kxibk6b5%2bvbjv0ikmzrxseb621wgp6l%2f4aoakg%2fnkxcgah3t6omkqpxw5dxjw4v27zcxpcu52gkqefbrhhg%3d%3d&x-amz-algorithm=aws4-hmac-sha256&x-amz-date=20210124t022835z&x-amz-signedheaders=host&x-amz-expires=14400&x-amz-credential=asiawsmx3snwskq7z4xt%2f20210124%2fus-west-2%2fs3%2faws4_request&x-amz-signature=c683f06894332a60da509cf935ef198efa0308e35cbb7f99f9be498f9daeda32"
  
  email = readline(
    "please provide your email (used to comply with terms of service for the nominatim geocoder):\n"
  )
  
  downloadDryadData = readline(
    "do you want to download a fresh copy of the data? Y = download the data directly from datadryad.org, N = enter your own file path:\n"
  )
  
  if ("Y" == downloadDryadData || "y" == downloadDryadData) {
    printMessageWithTimestamp("downloading from datadryad.org is not implemented at this time. exiting.")
    
    # download.file(url, destfile)
  } else {
    filepath = readline(
      "please provide the fully-qualified path to your local csv file. Leave blank to use the smaller sample data set.\n"
    )
    
    if ("" == filepath) {
      filepath = "data-raw/geocode_data/sample_zika_data.csv"
    }
    
    outputPath = paste0(tools::file_path_sans_ext(filepath), "_annotated.csv")
    
    printMessageWithTimestamp(paste0("Processing data in ", filepath))
    printMessageWithTimestamp(paste0("Writing results to ", outputPath))
    
    processFile(filepath, outputPath, email)
  }
}

executeInteractiveScript()
