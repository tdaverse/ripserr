# save input for 2dim test as RDS (after converting format, etc.)
test_input <- read.csv("data-raw/test_input_2dim.csv")

mat_input <- test_input$value
dim(mat_input) <- c(max(test_input$dim1),
                    max(test_input$dim2))

saveRDS(mat_input, file = "tests/testthat/input_2dim.rds")

# save output for 2dim test as RDS
