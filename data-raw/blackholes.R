library(magick)

powehi_img  <- "https://upload.wikimedia.org/wikipedia/commons/thumb/4/4f/Black_hole_-_Messier_87_crop_max_res.jpg/250px-Black_hole_-_Messier_87_crop_max_res.jpg"
sagA_img <- "https://upload.wikimedia.org/wikipedia/commons/thumb/9/96/EHT_Saggitarius_A_black_hole.tif/lossy-page1-240px-EHT_Saggitarius_A_black_hole.tif.jpg"
# read image and convert it to grayscale
powehi_img <- image_read(powehi_img) |> image_convert(colorspace = "gray")
sagA_img <- image_read(sagA_img) |> image_convert(colorspace = "gray")
# make the image into a bitmap matrix
powehi_mat <- as.numeric(image_data(powehi_img)) |> drop() |> t()
powehi <- powehi_mat[, ncol(powehi_mat):1]

sagA_mat <- as.numeric(image_data(sagA_img)) |> drop() |> t()
sagA <- sagA_mat[, ncol(sagA_mat):1]

#image(powehi, col = gray.colors(256))
#image(sagA, col = gray.colors(256))

# as.numeric() needed to change from raw datatype
# t() and vertical flip needed as {Magick} works from bottom row up, but 
#`image()` plots where the x-axis is row index and the y-axis is the column index
# see ?image

save(powehi, file = "powehi.rda")
save(sagA, file = "sagA.rda")
