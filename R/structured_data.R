#' @title \emph{structured_data}: function to organize structured datasets.
#' @description The function is used to convert a collection of structured datasets into an object which can be used in \code{speciesModel}. The three types of structured data allowed by this function are present absence (PA), present only (PO) and counts/abundance datasets, which are controlled using the \code{datasetType} argument. The other arguments of this function are used to specify the appropriate variable (such as response name, trial name, species name and coordinate name) names in these datasets.
#'
#' @param ... The datasets used in the model. May be either datasets or a SpatialPointsDataFrames.
#' @param datasetType A vector which gives the type of dataset. Must be either \code{'count'}, \code{'PO'} or \code{'PA'}. Defaults to \code{NULL}.
#' @param responsePA Name of the response variable in the \code{PA} datasets. Defaults to \code{NULL}.
#' @param trialsPA Name of the trial name variable in the \code{PA} datasets. Defaults to \code{NULL}.
#' @param responseCount Name of the response variable in the \code{count} datasets. Defaults to \code{NULL}.
#' @param speciesName Name of the species variable name. Defaults to \code{'species'}.
#' @param coordinateNames Name of the coordinates used in the model. Defaults to \code{c('longitude', 'latitude')}.
#'
#' @import methods
#' @import sp
#'
#' @examples {
#'
#'  #Objects required for example
#'  data("PA_redlist")
#'  speciesNames <- c('Fraxinus excelsior', 'Ulmus glabra')
#'
#'  #Set up structured dataset
#'    dataObj <- structured_data(PA_redlist, datasetType = c('PA'), responsePA = 'individualCount',
#'                              speciesName = 'species',
#'                              coordinateNames = c("longitude", "latitude" ))
#'
#' }
#'
#' @return An S4 object of class \code{structuredData}, which contains three slots for the data of each observation model allowed in the framework.
#' @export

structured_data <- function(..., datasetType = NULL, responsePA = NULL,
                          trialsPA = NULL, responseCount = NULL,
                          speciesName = 'species',
                          coordinateNames = c('longitude', 'latitude')) {

  if (is.null(datasetType)) stop('datasetType cannot be NULL.')

  data <- list(...)

  if (!all(datasetType%in%c('count', 'PO','PA'))) stop('datasetType must be a vector with values: "count", PO" or "PA".')

  if (length(data) != length(datasetType)) stop('The length of the datasets is not equal to the length of datasetType.')

  if ('PA'%in%datasetType && is.null(responsePA)) stop('PA datasets are included in datasetType, but responsePA is NULL.')

  if ('count'%in%datasetType && is.null(responseCount)) stop ('Count datasets are included in datasetType but responseCount is NULL.')

  if (is.null(speciesName)) stop('Please supply a species variable name.')

  if (is.null(coordinateNames) | length(coordinateNames) != 2) stop('coordinateNames has to be a vector of two elements describing the coordinates used in the model.')

  checkCoords <- sapply(data, function(x) {

    if (inherits(x, 'data.frame')) {

      all(coordinateNames%in%names(x))

    }
    else
      if (inherits(x,'Spatial')) {

        all(coordinateNames%in%colnames(x@coords))

      }


  })

  if (!all(checkCoords)) stop('All datasets are required to have the same coordinate names specified by coordinateNames.')

  checkSpecies <- sapply(data, function(x) {

    if (inherits(x, 'data.frame')) {

      speciesName%in%names(x)

    }
    else
      if (inherits(x,'Spatial')) {

     speciesName%in%names(x@data)

      }


  })

  if (!all(checkSpecies)) stop('All datasets names are required to have the same species variable name specified with speciesName.')

  data <- lapply(data, function(data) {

    if (inherits(data, 'data.frame')) {

      data[,speciesName] <- gsub(' ','_', data[,speciesName])

      return(data)
    }

    else {

      data@data[,speciesName] <- gsub(' ', '_', data@data[,speciesName])

      return(data)

    }

  })


  data <- sapply(data, function(data) {

  if (inherits(data, 'data.type')) {

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

  PO_data <- data[datasetType == 'PO']
  PA_data <- data[datasetType == 'PA']
  Count_data <- data[datasetType == 'count']

  if (length(PO_data) == 0) PO_data <- NULL
  if (length(PA_data) == 0) PA_data <- NULL
  if (length(Count_data) == 0) Count_data <- NULL

  data_object <- new('structuredData',
                     dataPO = PO_data,
                     dataPA = PA_data,
                     dataCount = Count_data)


  attr(data_object, 'datasetType') <- datasetType

  if (is.null(responsePA)) responsePA <- 'responsePA'

  attr(data_object, 'responsePA') <- responsePA

  if (is.null(trialsPA)) trialsPA <- 'trialsPA'

  if (is.null(responseCount)) responseCount <- 'responseCount'

  attr(data_object, 'responseCount') <- responseCount

  attr(data_object, 'trialsPA') <- trialsPA
  attr(data_object, 'speciesName') <- speciesName
  attr(data_object, 'coordinateNames') <- coordinateNames

  data_object

}
