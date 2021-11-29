#' Function to convert a collection of datasets into an object which can be used in \code{speciesModel}.
#'
#' @param ... The datasets used in the model. May be either datasets or a SpatialPointsDataFrames.
#' @param datasetType A vector which gives the type of dataset. Must be either \code{'PO'} or \code{'PA'}. Defaults to \code{NULL}.
#' @param responsePA Name of the response variable in the \code{PA} datasets. Defaults to \code{NULL}.
#' @param trialsPA Name of the trial name variable in the \code{PA} datasets. Defaults to \code{NULL}.
#' @param speciesName Name of the species variable name. Defaults to \code{'species'}.
#' @param coordinateNames Name of the coordinates used in the model. Defaults to \code{c("latitude","longitude")}.
#'
#' @export

structured_data <- function(..., datasetType = NULL, responsePA = NULL,
                          trialsPA = NULL, speciesName = 'species',
                          coordinateNames = c('latitude', 'longitude')) {

  if (is.null(datasetType)) stop('datasetType cannot be NULL.')

  data <- list(...)

  if (!all(datasetType%in%c('PO','PA'))) stop('datasetType must be a vector with values: "PO" or "PA".')

  if (length(data) != length(datasetType)) stop('The length of the datasets is not equal to the length of datasetType.')

  if ('PA'%in%datasetType & is.null(responsePA)) stop('PA datasets are included in datasetType, but responsePA is NULL.')

  if (is.null(speciesName)) stop('Please supply a species variable name.')

  if (is.null(coordinateNames) | length(coordinateNames) != 2) stop('coordinateNames has to be a vector of two elements describing the coordinates used in the model.')

  checkCoords <- sapply(data, function(x) {

    if (class(x) == 'data.frame') {

      all(coordinateNames%in%names(x))

    }
    else
      if (inherits(x,'Spatial')) {

        all(coordinateNames%in%colnames(x@coords))

      }


  })

  if (!all(checkCoords)) stop('All datasets are required to have the same coordinate names specified by coordinateNames.')

  checkSpecies <- sapply(data, function(x) {

    if (class(x) == 'data.frame') {

      speciesName%in%names(x)

    }
    else
      if (inherits(x,'Spatial')) {

     speciesName%in%names(x@data)

      }


  })

  if (!all(checkSpecies)) stop('All datasets names are required to have the same species variable name specified with speciesName.')


  data <- lapply(data, function(data) {

    if (class(data) == 'data.frame') {

      data[,speciesName] <- gsub(' ','_', data[,speciesName])

      return(data)
    }

    else {

      data@data[,speciesName] <- gsub(' ', '_', data@data[,speciesName])

      return(data)

    }

  })


  data <- sapply(data, function(data) {

  if (class(data) == 'data.type') {

    if (length(names(data)) > 2) { ## ie coordinate names ++

  convertedData <- sp::SpatialPointsDataFrame(coords = data[,coordinateNames],
                                              data = data.frame(data[,!coordinateNames]),
                                              proj4string = CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'))
  names(convertedData@data)[, names(convertedData@data) == speciesName] <- 'species'
  return(convertedData)

    }
  else {

    convertedData <- sp::SpatialPointsDataFrame(coords = data[,coordinateNames],
                                                data = data.frame(species = data[,!coordinateNames]),
                                                proj4string = CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'))
    return(convertedData)


  }


  }

  else
    if (inherits(data, 'Spatial')) {

    proj4string(data) <- CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
    names(data@coords) <- coordinateNames
    return(data)

    }


  })

  data_names <- setdiff(as.character(match.call(expand.dots = TRUE)),
                        as.character(match.call(expand.dots = FALSE)))

  names(data) <- data_names

  if (any(datasetType == 'PO') & any(datasetType == 'PA')) {

  PO_data <- data[datasetType == 'PO']
  PA_data <- data[datasetType == 'PA']

  }
  else
    if (all(datasetType == 'PO')) {

    PO_data <- data
    PA_data <- NULL

    }
  else {

    PO_data <- NULL
    PA_data <- data

  }

  data_object <- new('structuredData',
                     dataPO = PO_data,
                     dataPA = PA_data)


  attr(data_object, 'datasetType') <- datasetType

  if (is.null(responsePA)) responsePA <- 'responsePA'

  attr(data_object, 'responsePA') <- responsePA

  if (is.null(trialsPA)) trialsPA <- 'trialsPA'

  attr(data_object, 'trialsPA') <- trialsPA
  attr(data_object, 'speciesName') <- speciesName
  attr(data_object, 'coordinateNames') <- coordinateNames

  return(data_object)

}
