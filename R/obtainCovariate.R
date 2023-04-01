#' @description
#' Function to obtain covariate layers for a study area
#' @param covariates A vector of covariate names to obtain.
#' @param countries A vector of country names to obtain covariates for.
#' @param projection Coordinate reference system to use in analysis.

obtainCovariate <- function(covariates, countries,
                            projection, path) {

  covariateLayers <- geodata::worldclim_country(country = countries,
                                                var = covariates,
                                                path = path)

  covariateLayers <- terra::project(covariateLayers, projection)

  covariateLayers

  }
