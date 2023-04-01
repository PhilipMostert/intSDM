#' @title \emph{species_model}: function to construct an integrated species distribution model, as well as other useful outputs from the model.
#' @description This function is used to construct integrated species distribution models using data specified with \code{\link{structured_data}} as well as present only data obtained from the Global Biodiversity Information Facility (GBIF).
#'
#' @param speciesNames A vector of species' names to collect from GBIF.
#' @param date Vector of length two denoting the date range to select species from. Defaults to \code{NULL}.
#' @param gbifOpts A named list of additional options to filter the GBIF records. See \code{?rgbif::occ_search} for the filters which may be applied to the GBIF records. Defaults to \code{list(coordinateUncertaintyInMeters = '0,1000')} which specifies a desired range of coordinates uncertainty (in meters) between 0 and 1000.
#' @param structuredData Additional datasets to integrate with the presence only GBIF data. See the \code{structured_data} function. Defaults to \code{NULL}.
#' @param spatialCovariates Spatial covariates to include in the model. May be a \code{Raster} or \code{Spatial} object. Cannot be non-\code{NULL} if \code{worldclimCovariates} is non-\code{NULL}.
#' @param worldclimCovariates Names of the covariates to extract from Worldclim. Defaults to \code{NULL}; cannot be non-\code{NULL} if \code{spatialCovariates} is non-\code{NULL}.
#' @param res Resolution for the WorldClim covariates. Valid values are: \code{0.5,2.5,5,10}. Defaults to \code{0.5}.
#' @param scale Should the spatial covariates be scaled. Defaults to \code{FALSE}.
#' @param location Which area of Norway to model. Defaults to \code{'Norway'} which suggests a model for the entire county.
#' @param boundary SpatialPolygons object of the study area. If \code{NULL} an object may be formed with \code{location}.
#' @param return Object to return. Has to be one of \code{c('boundary', 'species', 'species plot', 'mesh', 'mesh plot', 'model', 'predictions', 'predictions map')}.
#' @param mesh An inla.mesh object to include in the model. Defaults to \code{NULL}.
#' @param meshParameters A list of inla.mesh arguments to create a mesh if \code{mesh = NULL}.
#' @param spdeModel \code{inla.spde} model used in the model. May be a named list where the name of the spde object is the name of the associated dataset. Default \code{NULL} uses \code{inla.spde2.matern}.
#' @param biasField Should a second (bias) random field be added to the GBIF data. Defaults to \code{FALSE}.
#' @param biasModel \code{inla.spde} model used for the bias field. Requires \code{biasModel} to be \code{TRUE}. Default \code{NULL} uses \code{inla.spde2.matern}.
#' @param projection CRS projection to use. Defaults to \code{CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')}.
#' @param limit Set the number of species downloaded. Defaults to \code{10000}.
#' @param options A list of \code{INLA} and \code{inlabru} options. Defaults to \code{NULL}.
#' @param ... Additional arguments used in PointedSDMs's intModel function.
#'
#' @import ggplot2
#' @import sp
#' @import methods
#' @import stats
#' @import inlabru
#' @import PointedSDMs
#' @importFrom raster getData
#' @importFrom raster crop
#' @importFrom raster merge
#' @importFrom raster mask
#@importFrom maps map
#' @importFrom spocc occ
#' @importFrom dplyr bind_rows
#@importFrom maptools map2SpatialPolygons
#'
#' @examples {
#'
#'
#' \dontrun{
#'  if(requireNamespace('INLA')) {
#'
#'  #Objects required for example
#'  data("PA_redlist")
#'  speciesNames <- c('Fraxinus excelsior', 'Ulmus glabra')
#'
#'  #Set up structured dataset
#'    dataObj <- structured_data(PA_redlist, datasetType = c('PA'), responsePA = 'individualCount',
#'    speciesName = 'species',
#'    coordinateNames = c("longitude", "latitude" ))
#'
#'   #Get species map
#'     predictions <- species_model(return = 'species plot',
#'     boundary = boundary, speciesNames = species,
#'     limit = 10, structuredData = dataObj,
#'     meshParameters = list(cutoff=0.08, max.edge=c(1, 3), offset=c(1,1)),
#'     worldclimCovariates ='Annual Mean Temperature')
#'   }
#'
#'   }
#'
#' }
#'
#' @return The return of the function is determined by the argument \code{return}. For the different values of return: \item{\code{boundary}}{A \code{SpatialPolygon} of the boundary used in the model,} \item{\code{species}}{A \code{data.frame} object of the species' coordinates used in the model,} \item{\code{species plot}}{a \code{ggplot} plot of the species across a map,} \item{\code{mesh}}{an \code{inla.mesh} object,} \item{\code{mesh plot}}{a plot of the \code{inla.mesh object},} \item{model}{the integrated model,} \item{\code{predictions}}{predictions from the integrated model (on the linear scale),} \item{\code{predictions map}}{A \code{ggplot} plot of the predictions across the boundary.}
#'
#' @export

species_model <- function(speciesNames,
                          date = NULL,
                          gbifOpts = list(coordinateUncertaintyInMeters = '0,1000'),
                          structuredData = NULL,
                          spatialCovariates = NULL,
                          worldclimCovariates = NULL,
                          res = 0.5, scale = FALSE,
                          location = 'Norway', boundary = NULL,
                          return = 'predictions map', mesh = NULL,
                          meshParameters = NULL, spdeModel = NULL,
                          biasField = FALSE, biasModel = NULL,
                          projection = CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'),
                          limit = 10000, options = list(), ...) {

  if (missing(speciesNames) && is.null(structuredData) && return %in% c('model', 'predictions', 'predictions map')) stop('Please provide at least one of speciesNames and structuredData.')

  if (!is.null(structuredData)) {

    if (!inherits(structuredData, 'structuredData')) stop('Please run your additional data through the structured_data function.')

  }

  worldclimIndex <- data.frame(index = paste0('layer.',1:19),
                               variable = c('Annual Mean Temperature', 'Mean Diurnal Range',
                                            'Isothermality', 'Temperature Seasonality', 'Max Temperature of Warmest Month',
                                            'Min Temperature of Coldest Month', 'Temperature Annual Range',
                                            'Mean Temperature of Wettest Quarter', 'Mean Temperature of Driest Quarter',
                                            'Mean Temperature of Warmest Quarter', 'Mean Temperature of Coldest Quarter',
                                            'Annual Precipitation', 'Precipitation of Wettest Month',
                                            'Precipitation of Driest Month', 'Precipitation Seasonality ',
                                            'Precipitation of Wettest Quarter', 'Precipitation of Driest Quarter',
                                            'Precipitation of Warmest Quarter', 'Precipitation of Coldest Quarter'))

  if (length(return) > 1) stop('return must contain only one element.')

  if (!return%in%c('boundary','species', 'species plot', 'mesh', 'mesh plot', 'model', 'predictions', 'predictions map')) stop('return is not one of: boundary, species, species plot, mesh, mesh plot, model, predictions, predictions map.')

  #if (!is.null(boundary)) {

  #  if(class(boundary) != 'SpatialPolygons' | class(boundary) != 'SpatialPolygonsDataFrame') stop('Boundary needs to be a SpatialPolygon*s object.')

  #}

  if (!is.null(worldclimCovariates) & !is.null(spatialCovariates)) stop("One of worldclimCovariates and spatialCovariates must be NULL. Please only chose one set of spatial covariates to model (or none)")

  if (!is.null(spatialCovariates)) {

    if (inherits(spatialCovariates, 'data.frame')) stop('data.frame spatial covariates are not permitted. Please convert these to a Raster* or Spatial* object.')
  }

  if (!is.null(worldclimCovariates)) {

    if (!any(worldclimCovariates%in%worldclimIndex[,2])) stop('worldclimCovariates given includes variable names not present in the worldclim dataset.')

  }

  if (is.null(boundary)) {
    ## or get from getData('GADM' , country="NOR", level=1) using RASTER ??
    ## probably better but super slow
    #boundary <- fhimaps::norway_nuts3_map_b2020_default_sf
    #boundary <- as(boundary,'Spatial')
     if (location %in% c('Norway', 'norway')) {

        if (file.exists('~/Data-raw')) norwayfill <- geodata::gadm(country = 'Norway', path = '~/Data-raw', resolution = 2)
        else norwayfill <- geodata::gadm(country = 'Norway', path = getwd(), resolution = 2)

        boundary <- as(norwayfill, 'Spatial')
        proj4string(boundary) <- projection

       #norwayfill <- maps::map("world", "norway", fill=TRUE, plot=FALSE,
       #                   ylim=c(58,72), xlim=c(4,32))
       #IDs <- sapply(strsplit(norwayfill$names, ":"), function(x) x[1])
       #boundary <- maptools::map2SpatialPolygons(norwayfill, IDs = IDs,
       #                                   proj4string = projection)

    }
    else {

      #stop('TODO::Subset to get the separate counties through mapNO[[1]][1:11]. Make data.frame with index for all counties.')

      boundary <- raster::getData(name = 'GADM', country = 'NOR', level = 1)
#do paste('\uxxx)
      if (!any(location%in%c("Akershus", paste0('\u00C3',"stfold"), "Aust-Agder", "Buskerud", "Finnmark", "Hedmark", "Hordaland",
                              paste0("M", '\u00F8' ,"re og Romsdal"), paste0("Nord-Tr", '\u00F8',"ndelag"), "Nordland", "Oppland", "Oslo", "Rogaland", "Sogn og Fjordane",
                              paste0("S", '\u00F8',"r-Tr", '\u00F8',"ndelag"), "Telemark", "Troms", "Vest-Agder", "Vestfold"))) stop(paste('At least one of the locations provided is not a valid county in Norway. NOTE:', paste0('Tr', '\u00F8', 'ndelag'), 'is given as', paste0('Nord-Tr', '\u00F8', 'ndelag'), 'and', paste0("S", '\u00F8',"r-Tr", '\u00F8',"ndelag")))

      warning('Location is given as a region of Norway. Mesh creation may be slow.')

      boundary <- boundary[boundary$NAME_1 == location,]
    }

    if (return == 'boundary') {

      return(boundary)

    }

  }

  if (!is.null(structuredData)) {
    ##Rename all coord names to "latitude";"longitude" to reduce arguments required
    ##Also standardize speciesName to 'species'

    coordinateNames <- unlist(attributes(structuredData)['coordinateNames'])
    speciesName <- unlist(attributes(structuredData)['speciesName'])
    responsePA <- unlist(attributes(structuredData)['responsePA'])
    trialsPA <- unlist(attributes(structuredData)['trialsPA'])
    responseCount <- unlist(attributes(structuredData)['responseCount'])

    structuredData <- append(structuredData@dataPO, append(structuredData@dataPA, structuredData@dataCount))

    if (any(coordinateNames != c('longitude', 'latitude'))) {

      structuredData <- lapply(structuredData, function(x) {

        colnames(x@coords) <- c('longitude', 'latitude')
        return(x)

      })

      if (speciesName != 'species') {

        structuredData <- lapply(structuredData, function(x) {

          names(x@data)[names(x) == speciesName] <- 'species'
          return(x)

        })

        #responsePA <- attributes(structuredData)['responsePA']

        if (!is.null(responsePA)) responsePA <- 'responsePA'

        #trialsPA <- attributes(structuredData)['trialsPA']

        if (!is.null(trialsPA)) trialsPA <- 'trialsPA'

      }

    }
  }
  else {
    speciesName <- 'species'
    responsePA <- 'responsePA'
    trialsPA <- 'trialsPA'
    responseCount <- 'responseCount'

  }

  if (is.null(mesh)) {
    if (is.null(meshParameters)) stop('meshParameters cannot be NULL if mesh is NULL.')
     if (any(!c("cutoff", "max.edge", "offset") %in% names(meshParameters))) stop("'cutoff', 'max.edge' and 'offset' need to be in meshParameters")
    boundary <- as(boundary, 'SpatialPolygons')
    message('Making inla.mesh object:')
    mesh <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(boundary),
                           cutoff = meshParameters$cutoff,
                           max.edge = meshParameters$max.edge,
                           offset = meshParameters$offset)
    mesh$proj4string <- projection
    mesh$crs <- projection

    if (return == 'mesh') {

      return(mesh)

    }

  }
  else {
    ###TODO
    ##Make BBOX from mesh object

  }

  if (return == 'mesh plot') {

    meshPlot <- ggplot() +
      gg(mesh) +
      ggtitle('Plot of mesh') +
      xlab('Longitude') +
      ylab('Latitude') +
      theme(plot.title = element_text(hjust = 0.5))
    return(meshPlot)

  }
  ## Should this be an optional argument???
  if (!missing(speciesNames)){

  message('Obtaining GBIF species data:')
##Mighht be worth doing this as a R6 + adding argument related to occuranceStatus and individualCount (note NA)
  if (is.null(date)) {

    species_data <- spocc::occ(query = speciesNames,
                               limit = limit,
                               gbifopts = gbifOpts,
                               geometry = boundary@bbox,
                               has_coords = TRUE)

  }
  else  species_data <- spocc::occ(query = speciesNames,
                             limit = limit,
                             gbifopts = gbifOpts,
                             geometry = boundary@bbox,
                             has_coords = TRUE,
                             date = date)

  species_data <- do.call(dplyr::bind_rows,
                          species_data$gbif$data)
  ##Why am I getting species other than the ones from speciesNames????
  ##Or think of nice ways to standardize all species names::
  ##Somehow even relate it back to structuredData
  ##But is tough since species in structuredData does not
  ##Necessarily need to be in the PO data and vice versa ...
  species_data$species <- gsub(' ', '_', species_data$species)

  all_data <- sp::SpatialPointsDataFrame(coords = cbind(species_data$longitude, species_data$latitude),
                                         data = data.frame(species = species_data$species),
                                         proj4string = projection)
  ##Correct??
  colnames(all_data@coords) <- c('longitude', 'latitude')

  ##Select species only within boundary???
  ##Why are species leaking through the bbox??
  all_data <- all_data[boundary, ]
  all_data <- list(dataGBIF = all_data)

  ###Things to do::
  ## so structured data is a list of sp objects -- possibly causing issues
  ## Create a list; all data with one element; all_data
  ## Then for loop all elements from structured data::
  ## Note that names may not be ket as is?
  ######
  if (!is.null(structuredData)) {

    for (dataset in 2:(length(structuredData) + 1)) {

      all_data[[dataset]] <- structuredData[[(dataset - 1)]]
      names(all_data)[[dataset]] <- names(structuredData)[[(dataset - 1)]]

    }

  }

  }
  else all_data <- structuredData

  if (return == 'species') return(all_data)

  if (return == 'species plot') {

    data_to_plot <- lapply(all_data, function(x) x[, 'species'])
    data_to_plot <- do.call(sp::rbind.SpatialPointsDataFrame, data_to_plot)
    species <- 'species'
    plot <- ggplot() + gg(boundary)  + gg(data_to_plot, aes(colour = species))
    return(plot)
  }


  if (!is.null(worldclimCovariates)) {
    ##Some reason worldclim includes half of Norway in the one lat and the other half in the other lat when res = 0.5
    ##Need to glue the two raster files together
    message('Obtaining worldclim covariates:')
    bioclimS <- raster::getData("worldclim", var = "bio", res = res, lon = 5, lat = 60)
    bioclimN <- raster::getData("worldclim", var = "bio", res = res, lon = 5, lat = 70)

    r1 <- raster::crop(bioclimN, bbox(boundary))
    r2 <- raster::crop(bioclimS, bbox(boundary))

    mergedNO <- raster::merge(r1, r2)
    spatialCovariates <- raster::mask(mergedNO, boundary)

    #spatialCovariates <- subset(spatialCovariates,worldclimIndex[,'index'][worldclimIndex['variable'] == worldclimCovariates])
    spatialCovariates <- as(spatialCovariates, 'SpatialPixelsDataFrame')
    spatialCovariates <- spatialCovariates[,worldclimIndex[,'index'][worldclimIndex['variable'] == worldclimCovariates]]
    names(spatialCovariates) <- gsub(' ', '_', worldclimCovariates)

  }
  else
    if (is.null(spatialCovariates) & is.null(worldclimCovariates)) spatialCovariates <- NULL

  ##Don't know why this needs to be a dataframe...
  if (!is.null(spatialCovariates) & scale) spatialCovariates@data <- data.frame(scale(spatialCovariates@data))

  message('Organizing the data:')


  organized_data <- PointedSDMs::intModel(all_data, spatialCovariates = spatialCovariates, Coordinates = c('longitude', 'latitude'),
                                        Mesh = mesh, responseCounts = responseCount, responsePA = responsePA,
                                        trialsPA = trialsPA, speciesName = 'species',
                                        Projection = projection, ...)

  if (!is.null(spdeModel)) {

    if (!is.null(organized_data$.__enclos_env__$private$Spatial[['sharedField']])) organized_data$spatialFields$sharedField[['sharedField']] <- spdeModel
    else stop('spdeModel provided but no shared spatial effect included in the model.')
  }

  if (biasField) organized_model$addBias(datasetName = dataGBIF, biasField = biasModel)
  ##Print model components here? new argument called verbose or something?

  message('Running model:')

  #spatialModel <- PointedSDMs::bru_sdm(data = organized_data, spatialcovariates = spatialCovariates,
  #                                    sharedspatial = TRUE, specieseffects = TRUE, spdemodel = spdeModel,
  #                                    options = options)

  spatialModel <- PointedSDMs::fitISDM(organized_data, options = options)

  if (return == 'model') {

    return(spatialModel)

  }

  ##Then predict
  message('Predicting model:')

  #covariateNames <- ifelse(is.null(spatialCovariates), NULL, names(spatialCovariates))
  mesh$crs <- projection
  ##Try calling predict.bru_sdm manually???
  ##Need to add something -- if # species == 1 then predict only spatial field and not species = TRUE?
  modelPredict <- predict(spatialModel, mesh = mesh, mask = boundary,
                                               spatial = TRUE, intercepts = TRUE,
                                               covariates = spatialModel[['spatCovs']][['name']],
                                               fun = '')

  if (return == 'predictions') {

    return(modelPredict)

  }

  message('Plotting predictions')
  predictPlot <- plot(modelPredict, whattoplot = 'mean',
                      plot = FALSE)

  predictPlot

}

