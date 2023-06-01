#' @description
#' Function to obtain covariate layers from _WorldClim_ around a specified area.
#' @param covariates A vector of covariate names to obtain.
#' @param countries A vector of country names to obtain covariates for.
#' @param projection Coordinate reference system to use in analysis.
#'
#' @return A \code{spatialRaster} object of the covariates across the specified area.

obtainCovariate <- function(covariates, countries,
                            projection, path) {

  ##Do for other covariate layers

  covariateLayers <- geodata::worldclim_country(country = countries,
                                                var = covariates,
                                                path = path)

  covariateLayers <- terra::project(covariateLayers, projection)

  covariateLayers

  }
