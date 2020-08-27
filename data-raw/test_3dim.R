# read and format input file
con_fin <- file("data-raw/test_input_3dim.txt", "r")
temp_input <- as.numeric(readLines(con = con_fin))
close(con_fin)

dim <- temp_input[1]
lens <- temp_input[2 : (1 + dim)]

test_input <- temp_input[-(1:(1 + dim))]
dim(test_input) <- lens

# save input file
saveRDS(test_input, "tests/testthat/input_3dim.rds")

# convert output from CSV to RDS
test_output <- read.csv("data-raw/test_output_3dim.csv")
saveRDS(test_output, "tests/testthat/output_3dim.rds")