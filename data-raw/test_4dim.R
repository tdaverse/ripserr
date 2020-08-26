# convert input data from CSV to RDS
input_data <- read.csv("raw-data/test_input_4dim.csv")
vals <- numeric(nrow(input_data))
dim(vals) <- c(max(input_data$dim1),
               max(input_data$dim2),
               max(input_data$dim3),
               max(input_data$dim4))

for (i in seq_len(nrow(input_data))) {
  vals[input_data$dim1[i], input_data$dim2[i], input_data$dim3[i], input_data$dim4[i]] <- input_data$value[i]
}

saveRDS(vals, file = "tests/testthat/input_4dim.rds")

# convert output data from CSV to RDS
output_data <- read.csv("raw-data/test_output_4dim.csv")
saveRDS(output_data, file = "tests/testthat/output_4dim.rds")
