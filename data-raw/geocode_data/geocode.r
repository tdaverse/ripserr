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

nominatimSleepDelay = 1.1

reverseGeocode <- function (longitude, latitude) 
{
  # Pause for 1.1 seconds, to honor the request of the Nominatim API documentation to not make more than one request per second.
  Sys.sleep(nominatimSleepDelay)
  
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
  
  if (file.exists(outputFilePath)) {
    # Remove past run results
    print("Warning: deleting previous run results.")
    file.remove(outputFilePath)
  }
  
  header <- read.csv(inputFilePath, nrows=1, header=FALSE)
  lineCount = countLines(inputFilePath)
  estimatedProcessingMinutes = lineCount * nominatimSleepDelay / 60
  print(paste0("Based on your file size, this process will take approximately ", round(estimatedProcessingMinutes, 2), " minutes"))
  
  for (rowIndex in 1:lineCount) {
    outputRow <- read.csv(inputFilePath, skip=rowIndex, nrows=1, header=FALSE, col.names=header)
    
    x  <- outputRow$X
    y  <- outputRow$Y
    
    if (identical(x, logical(0)) || identical(y, logical(0))) {
      next
    }
    
    if (0 == rowIndex %% 10) {
      print(paste0("   Progress: ", rowIndex, " rows of ", lineCount, " have been processed."))
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


defaultdownloadurl = "https://uc3-s3mrt1001-prd.s3.us-west-2.amazonaws.com/3877a3a1-b8e2-4c47-9a55-a03572325842/data?response-content-disposition=attachment%3b%20filename%3ddoi_10.5061_dryad.47v3c__v2.zip&response-content-type=application%2fzip&x-amz-security-token=iqojb3jpz2lux2vjeakacxvzlxdlc3qtmijgmeqciblms%2biaicmw%2bvtf3t1ywxfkrewykec%2b6csg2qd7rfq4aibssye%2fqdipxrayzrbfrfzaaeofoiovyiclxhfldvmoncq9awji%2f%2f%2f%2f%2f%2f%2f%2f%2f%2f8beaaaddq1mtgynjkxnde1nyimm3ggolpzburlhpuckpedmwo5kddo1adxqnytewmfuz2svn0redp927zd30da2ec3zba9avwi33%2fwxr2arvd4vwn%2be%2bibfv9o6hhmpzvl91us9o2j6txctic1ok%2fqax%2bnzu1nc1tgmsxr6rqj5aaerpqizlx5my1fqclfldw%2bbi035meqhcjmnv%2fykklr0itxx%2fsan0jqi55vl2u2lh1b%2fyyplx5rd5hdi1%2brxmgqjckuhelacdt2yshhsteynz%2buqoox6wtdz2qh%2f%2fvk3vsorll%2blbfz9ktncy0vvhalr4yvhdy8pvyvwrz415idy7qp0wppkioq1k%2bj5p9qc8stfo53sivzlzxxdofdy8ir84y0qz9dfel4avimm2oqht8aaubhvfjkzsh0mfefqkp%2bvdfebj1wmrkkio4ng8lczrwbngnksxzewu1fopujnqdhashzpwtftan9vy7v494fchqu9g%2flj94ndvpnkxa8ls1gqhct1tawhzbuilwm9abwjedeekprswtpnbnwdhnf61kj3ysrfkioozuqvrimiqkwjjgzgay67ahcocftdmsye3fdn9fixittzp6v%2bydy4bzju63gncyrlei1vwr%2b5ispdemztaudzge3lxsse%2bclaeic3izgn6qbn16fplgi%2bfwb9lold%2bglcecyh0wice7mqqns5lwkgbnqvztciaphpkkzxumpemhzkduz97mgds5yip%2f2b2yt3xklvuedju9b%2bka24n2q%2bwihtsrrbd8cq738ev6op12wkygmbe7r7wmhywpgcm1kxibk6b5%2bvbjv0ikmzrxseb621wgp6l%2f4aoakg%2fnkxcgah3t6omkqpxw5dxjw4v27zcxpcu52gkqefbrhhg%3d%3d&x-amz-algorithm=aws4-hmac-sha256&x-amz-date=20210124t022835z&x-amz-signedheaders=host&x-amz-expires=14400&x-amz-credential=asiawsmx3snwskq7z4xt%2f20210124%2fus-west-2%2fs3%2faws4_request&x-amz-signature=c683f06894332a60da509cf935ef198efa0308e35cbb7f99f9be498f9daeda32"

email = readline("please provide your email (used to comply with terms of service for the nominatim geocoder):\n")
downloadDryadData = readline("do you want to download a fresh copy of the data? Y = download the data directly from datadryad.org, N = enter your own file path:\n")

print(email)
print(downloadDryadData)

if ("Y" == downloadDryadData || "y" == downloadDryadData) {
  print("downloading from datadryad.org is not implemented at this time. exiting.")

  # download.file(url, destfile)
} else {
  filepath = readline("please provide the fully-qualified path to your local csv file. Leave blank to use the smaller sample data set.\n")

  if ("" == filepath) {
    filepath = "data-raw/geocode_data/sample_zika_data.csv"
  }
  
  outputPath = paste0(tools::file_path_sans_ext(filepath), "_annotated.csv")
  
  print(paste0("Processing data in ", filepath))
  print(paste0("Writing results to ", outputPath))
  
  processFile(filepath, outputPath)
}
