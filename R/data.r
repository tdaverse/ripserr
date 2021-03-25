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
#' @source \url{http://dx.doi.org/10.5061/dryad.47v3c}
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
