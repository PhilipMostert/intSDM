#' @title obtainCovariate: Function to obtain a covariate later for a specified country.
#' @description
#' Function to obtain covariate layers from _WorldClim_ around a specified area.
#' @param covariates A vector of covariate names to obtain.
#' @param countries A vector of country names to obtain covariates for.
#' @param projection Coordinate reference system to use in analysis.
#' @param path The path where the covariate will be saved.
#'
#' @import geodata
#' @import terra
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
