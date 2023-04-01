#' @description
#' Function to obtain a \code{ADD} object around a specified area.
#' @param name A vector of names of countries used in the analysis.
#' @param projection Coordinate reference system used.
#'
#'
obtainArea <- function(names, projection) {

  world <- gisco_get_countries()

  if (!all(names %in% world$NAME_ENGL)) stop('At least one name provided not a valid country.')

  countryMap <- world[world$NAME_ENGL %in% names,]
  countryMap <- st_transform(countryMap, crs = projection)

  ##Maybe some warning here if countries provided are not next to each other geographically
   #What scale of areas? Could we go even finer or should that be specified by the user

  countryMap

}
