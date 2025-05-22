#' @title _Aedes aegypti_ occurrences in Brazil in 2013
#'
#' @description A geographic dataset of known occurrences of _Aedes aegypti_
#'   mosquitoes in Brazil, derived from peer-reviewed and unpublished literature
#'   and reverse-geocoded to states.
#'
#' @format A [tibble][tibble::tibble] of 4411 observations and 13 variables:

#' \describe{
#'   \item{vector}{species identification (_aegypti_ versus _albopictus_)}
#'   \item{occurrence_id}{unique occurrence identifier}
#'   \item{source_type}{published versus unpublished, with reference identifier}
#'   \item{location_type}{point or polygon location}
#'   \item{polygon_admin}{admin level or polygon size; -999 for point locations}
#'   \item{y}{latitudinal coordinate of point or polygon centroid}
#'   \item{x}{longitudinal coordinate of point or polygon centroid}
#'   \item{status}{established versus transient population}
#'   \item{state_name}{name of reverse-geolocated state}
#'   \item{state_code}{two-letter state code}
#' }

#'   
#' @source \doi{10.5061/dryad.47v3c}
#' @examples
#' 
#' # calculate persistence data for occurrences in Acre
#' acre_coord <- aegypti[aegypti$state_code == "AC", c("x", "y"), drop = FALSE]
#' acre_rips <- vietoris_rips(acre_coord)
#' plot.new()
#' plot.window(
#'   xlim = c(0, max(acre_rips$death)),
#'   ylim = c(0, max(acre_rips$death)),
#'   asp = 1
#' )
#' axis(1L)
#' axis(2L)
#' abline(a = 0, b = 1)
#' points(acre_rips[acre_rips$dim == 0L, c("birth", "death")], pch = 16L)
#' points(acre_rips[acre_rips$dim == 1L, c("birth", "death")], pch = 17L)
"aegypti"

#' @title State-level predictors of mosquito-borne illness in Brazil
#'
#' @description A data set of numbers of cases of Dengue in each state of Brazil
#'   in 2013 and three state-level variables used in a predictive model.
#'
#' @format A data frame of 27 observations and 4 variables:

#' \describe{
#'   \item{POP}{state population in 2013}
#'   \item{TEMP}{average temperature across state municipalities}
#'   \item{PRECIP}{average precipitation across state municipalities}
#'   \item{CASE}{number of state Dengue cases in 2013}
#' }

#'   
#' @source
#' \url{https://web.archive.org/web/20210209122713/https://www.gov.br/saude/pt-br/assuntos/boletins-epidemiologicos-1/por-assunto},
#' \url{http://www.ipeadata.gov.br/Default.aspx},
#' \url{https://ftp.ibge.gov.br/Estimativas_de_Populacao/},
#' `https://www.ibge.Goiasv.br/geociencias/organizacao-do-territorio/estrutura-territorial/15761-areas-dos-municipios.html?edicao=30133&t=acesso-ao-produto`
#'      
#' \describe{
#'  Data pre-processing:
#'  After acquiring data from above links, we converted any dataset 
#'  embedded in PDF format to CSV. Using carried functionalities in the CSV 
#'  file, we sorted all datasets alphabetically based on state names to make 
#'  later iterations more convenient. Also, we calculated the annual average 
#'  temperature and added to the original dataset where it was documented by 
#'  quarter.
#' }
"case_predictors"

#' @title Black hole in Messier 87 captured by the Event Horizon Telescope
#'
#' @description A data set containing a grayscale bitmap of the Messier 87 black 
#' hole with 250x250 resolution.
#' 
#' @format A 250x250 matrix containing pixels of evaluated between 0 and 1.
#' 
#' @source 
#' \urlhttps://commons.wikimedia.org/wiki/File:Black_hole_-_Messier_87_crop_max_res.jpg
#' 
#' \describe{
#' Data pre-processing:
#' After acquiring the 250x250 image from the above link we used the 
#' {[magick](https://cran.r-project.org/package=magick)} for processing. We 
#' first convert the image from RGB to grayscale using the default 
#' "perceptually-weighted" conversion. We then acquired the raw 3D array, 
#' converted the data type to numerical, and dropped the singleton channel 
#' dimension. We then transposed the matrix and vertically flipped it to align 
#' with how `image()` reads a matrix.
#' }
#' 
#' @examples
#' # compute the persistence homology of the Messier 87 black hole 
#' m87bh <- m87_black_hole
#' image(m87bh)
#' ph <- cubical(m87bh)
#' # based on the image, we expect one especially prominent persistent feature in 1D
#' ph$persistence <- ph$death - ph$birth
#' ph_1 <- ph[ph$dimension == 1, ]
#' ph_1 <- ph_1[order(-ph_1$persistence), ]
#' head(ph_1)
#' 
#' plot.new()
#' plot.window(
#'   xlim = c(0, max(ph$death)),
#'   ylim = c(0, max(ph$death)),
#'   asp = 1
#' )
#' axis(1L)
#' axis(2L)
#' abline(a = 0, b = 1)
#' points(ph[ph$dim == 0L, c("birth", "death")], pch = 16L)
#' points(ph[ph$dim == 1L, c("birth", "death")], pch = 17L)
"m87_black_hole"


