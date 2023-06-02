#' @title obtainArea: Function to obtain a boundry for a specified country.
#' @description
#' Function to obtain a \code{sf} boundary object around a specified country.
#' @param names A vector of names of countries used in the analysis.
#' @param projection Coordinate reference system used.
#' @param ... Additional arguments passed to \link[giscoR]{gisco_get_countries}.
#'
#' @import sf
#' @import giscoR
#'
#' @return An \code{sf} object of the boundary of the specified country.
#'
#'
obtainArea <- function(names, projection, ...) {

  world <- giscoR::gisco_get_countries(year = 2020, ...)

  if (!all(names %in% world$NAME_ENGL)) stop('At least one name provided not a valid country.')

  countryMap <- world[world$NAME_ENGL %in% names,]
  countryMap <- sf::st_transform(countryMap, crs = as.character(projection))

  ##Maybe some warning here if countries provided are not next to each other geographically
   #What scale of areas? Could we go even finer or should that be specified by the user

  countryMap

}
