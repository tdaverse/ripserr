library(magick)

url_img  <- "https://upload.wikimedia.org/wikipedia/commons/thumb/4/4f/Black_hole_-_Messier_87_crop_max_res.jpg/250px-Black_hole_-_Messier_87_crop_max_res.jpg"

# read image and convert it to grayscale
img <- image_read(url_img) |> image_convert(colorspace = "gray")

# make the image into a bitmap matrix
img_mat <- as.numeric(image_data(img)) |> drop() |> t()
img_mat <- img_mat[, ncol(img_mat):1]
# image(img_mat, col = gray.colors(256))

# as.numeric() needed to change from raw datatype
# t() and vertical flip needed as {Magick} works from bottom row up, but 
#`image()` plots where the x-axis is row index and the y-axis is the column index
# see ?image

save(img_mat, file = "Black_Hole-Messier_87.rda")
