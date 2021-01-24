library(RJSONIO)
library(hash)
library(R.utils)

brazilStateCodeLookup <- hash(
  "Acre"="AC",
  "Alagoas"="AL",
  "Amapá"="AP",
  "Amazonas"="AM",
  "Bahia"="BA",
  "Ceará"="CE",
  "Distrito Federal"="DF",
  "Espírito Santo"="ES",
  "Goiás"="GO",
  "Maranhão"="MA",
  "Mato Grosso"="MT",
  "Mato Grosso do Sul"="MS",
  "Minas Gerais"="MG",
  "Pará"="PA",
  "Paraíba"="PB",
  "Paraná"="PR",
  "Pernambuco"="PE",
  "Piauí"="PI",
  "Rio de Janeiro"="RJ",
  "Rio Grande do Norte"="RN",
  "Rio Grande do Sul"="RS",
  "Rondônia"="RO",
  "Roraima"="RR",
  "Santa Catarina"="SC",
  "São Paulo"="SP",
  "Sergipe"="SE",
  "Tocantins"="TO"
)

reverseGeocode <- function (longitude, latitude) 
{
  # Pause for 1.1 seconds, to honor the request of the Nominatim API documentation to not make more than one request per second.
  Sys.sleep(1.1)
  
  url <- paste0(
    "https://nominatim.openstreetmap.org/reverse?format=json&zoom=5&addressdetails=1&lon=",
    longitude,
    "&lat=", 
    latitude
    )
  
  response <- readLines(url, encoding = "UTF-8", warn = FALSE)
  response <- fromJSON(response)
  
  state <- tryCatch(response$address["state"], error = function(e) NULL)
  country <- tryCatch(response$address["country"], error = function(e) NULL)
  country_code <- tryCatch(response$address["country_code"], error = function(e) NULL)
  
  result <- list()
  result["state_name"] <- state
  result["state_code"] <- brazilStateCodeLookup[[state]]
  result["country_name"] <- country
  result["country_code"] <- country_code
  
  return(result)
}

processFile <- function(inputFilePath, outputFilePath)
{
  print("Reverse geocoding has started.")
  
  header <- read.csv(inputFilePath, nrows=1, header=FALSE)
  for (rowIndex in 1:countLines(inputFilePath)) {
    outputRow <- read.csv(inputFilePath, skip=rowIndex, nrows=1, header=FALSE, col.names=header)
    
    x  <- outputRow$X
    y  <- outputRow$Y
    
    if (identical(x, logical(0)) || identical(y, logical(0))) {
      next
    }
    
    if (0 == rowIndex %% 10) {
      print(paste0("Progress: row ", rowIndex, " of ", nrow(data), "."))
    }
    
    geocodeResult <- reverseGeocode(x, y)
    
    outputRow["GEOCODED_STATE_NAME"] <- geocodeResult$state_name
    outputRow["GEOCODED_STATE_CODE"] <- geocodeResult$state_code
    outputRow["GEOCODED_COUNTRY_NAME"] <- geocodeResult$country_name
    outputRow["GEOCODED_COUNTRY_CODE"] <- geocodeResult$country_code
    
    writeColumnNames <- (1 == rowIndex)
    
    try
    
    suppressWarnings(write.table(
      outputRow, 
      outputFilePath,
      append = TRUE,
      sep = ",",
      col.names = writeColumnNames,
      row.names = FALSE,
      quote = FALSE,
      fileEncoding="UTF8"))
  }
  
  print("Reverse geocoding has completed.")
}

# 
# defaultDownloadUrl = "https://uc3-s3mrt1001-prd.s3.us-west-2.amazonaws.com/3877a3a1-b8e2-4c47-9a55-a03572325842/data?response-content-disposition=attachment%3B%20filename%3Ddoi_10.5061_dryad.47v3c__v2.zip&response-content-type=application%2Fzip&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEAkaCXVzLXdlc3QtMiJGMEQCIBlmS%2BIAICmW%2BvTf3T1YWxfKrEwYkEc%2B6cSG2qd7RFQ4AiBSSYe%2FqdipxraYzrBFrFZaaEoFOIoVyIcLXhfLdVMonCq9Awji%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAAaDDQ1MTgyNjkxNDE1NyIMm3ggoLpzburLHPUCKpEDmWo5kdDO1aDXQNyTEWMFuZ2SVN0reDp927Zd30Da2ec3Zba9avwI33%2FWXr2ArvD4VWN%2Be%2BIbFv9O6Hhmpzvl91US9o2j6txctIc1OK%2FQaX%2Bnzu1NC1tGmSXr6Rqj5aaErPQIzLX5MY1FQClfldW%2BBi035mEqhCjMnV%2FyKKLR0ITXx%2FSaN0JQi55Vl2u2lH1B%2FyypLx5Rd5HDI1%2BrXmgqjCKUhELAcdT2YsHHSTeYnz%2BUQOOx6wTdZ2qh%2F%2Fvk3VSORLl%2BLBFz9kTncy0vVhalr4YvhDy8pvyVWRZ415idY7qp0wpPKiOq1k%2Bj5p9QC8stfo53sIvzLZXXdoFDy8iR84Y0qZ9dFEL4AViMm2oqHt8AaUBhVfJkZsH0mFefqKP%2BVdFeBj1WMrkkIo4Ng8lczrwBNGnkSxZewu1FoPUjnQdHashzPwtFTan9vY7V494FCHQU9g%2FLJ94NDvPnkxa8lS1GQHCt1taWhzBuiLwm9aBwJedEEkpRSwtPnBnwdhNf61kJ3YSRfKIOOzuQvrimiQkwjJGzgAY67AHCOcfTdMsYE3fdn9fIXITTZp6v%2BYdY4BzJU63gncyrlEi1Vwr%2B5ispDEmzTaudzgE3LxSSe%2BCLaeic3iZGN6qbN16fplGI%2BFwB9LOld%2BGlcEcYH0WICE7mqQNS5lWkGbnqVZtCiapHPkkzXUMPEmhZKduZ97mgdS5yiP%2F2B2Yt3XKLVueDjU9b%2BKA24n2q%2BwIHtsRrbd8CQ738Ev6op12WKYGmbe7r7wmhYWPGCM1KxiBk6b5%2BVbJV0iKmzRxseb621wGP6l%2F4AOAkg%2FnKxcgaH3t6OMKQPXW5dXJw4v27ZcXpcu52gKqEFBrHhg%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20210124T022835Z&X-Amz-SignedHeaders=host&X-Amz-Expires=14400&X-Amz-Credential=ASIAWSMX3SNWSKQ7Z4XT%2F20210124%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Signature=c683f06894332a60da509cf935ef198efa0308e35cbb7f99f9be498f9daeda32"
# 
# email = readline("Please provide your email (used to comply with TOS for the Nominatim geocoder):\n")  
# downloadDryadData = readline("Do you want to download a fresh copy of the data? Y = download the data directly from datadryad.org, N = enter your own file path:\n")
# 
# print(email)
# print(downloadDryadData)
# 
# if ("Y" == downloadDryadData || "y" == downloadDryadData) {
#   print("Downloading from datadryad.org is not implemented at this time. Exiting.")
#   
#   # download.file(url, destfile)
# } else {
#   filePath = readline("Please provide the fully-qualified path to your local csv file:\n") 
#   
#   print(filePath)
# }

processFile("./data-raw/geocode_data/sample_zika_data.csv", "./data-raw/geocode_data/sample_zika_data_annotated.csv")