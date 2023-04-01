##Redo species_model:

species_model <- R6::R6Class(classname = 'species_model', public = list(

  print = function(...) {

    cat('intSDM reproducible workflow summary:\n\n')

    cat('Research area: ')

    if (is.null(private$Area)) {

      cat('No research area provided: please add one using `.$addArea`.\n\n')

    } else {

      if (!is.null(private$Countries)) cat(private$Countries, '\n\n')
      else cat('Own boundry specified, research area unknown.\n\n')

      }

    cat('Datasets added: ')

    if (length(private$dataGBIF) == 0 && length(private$dataStructured) == 0) {

    cat('No species have been added to the workflow. Please add structured datasets using `.$addStructured` and data from GBIF using `.$addGBIF`.\n\n')


    } else {

      cat(names(c(private$dataGBIF, private$dataStructured)), '\n\n')



    }

    cat('Environmental covariates added: ')

    if (is.null(private$Covariates)) {

      cat('No covariates have been added to the workflow. Please add covariates using `.$addCovariates`.\n\n')

    } else {

      cat(names(private$Covariates), '\n\n')


    }

    cat('Model formula: ')

    if (is.null(private$CVMethod)) cat('No Cross-validation specified. Please specify using `.$crossValidation`.\n\n')
    else cat('Cross-validation:', paste(private$CVMethod, collapse = ', '),'\n\n')

    if (is.null(private$Output)) cat('No Output has been specified. Please specify using `.$workflowOutput`.')
    else cat('Model output:', paste(private$Output, collapse = ', '))

   }
  ,
  plot = function(Mesh = FALSE,
                  Boundary = TRUE,
                  Species = FALSE,
                  Covariates = FALSE) {

    if (sum(Mesh, Boundary, Species, Covariates) == 0) stop('At least one of Mesh, Boundary, Species or Covariates needs to be chosen.')

    if (Mesh) {

      if (is.null(private$Mesh)) stop('No mesh object provided. Please add one using `.$addMesh()`.')

      meshComponent <- gg(private$Mesh)

    } else meshComponent <- NULL

    if (Boundary) {

      if (is.null(private$Area)) stop('No boundary object proived. Please add one using `.$addArea()`.')

      boundaryComponent <- geom_sf(data = private$Area)

    } else boundaryComponent <- NULL

    if (Species) {

      if (all(is.null(private$dataGBIF), is.null(private$dataStructured))) stop('No data objects provided. Please add using either `.$addGBIF` or `.$addStructured`.')

      plotSpecies <- do.call(rbind, args = c(private$dataGBIF, private$dataStructured)) ## Need a variable with which to separate on ie species name

      speciesComponent <- geom_sf(data = plotSpecies) #+ facet_wrap( ~ speciesName) ##Need to convert this speciesName arg for all species



    } else speciesComponent <- NULL

    if (Covariates) {

      if (is.null(private$Covariates)) stop('No covariates object provided. Please add using `.$addCovariates()`.')

      plotList <- vector(mode = 'list',
                         length = length(private$Covariates))

      for (cov in names(private$Covariates)) {


        plotList[[cov]] <- ggplot() +
                             boundaryComponent +
                             meshComponent +
                             speciesComponent +
                             geom_spatraster(data = private$Covariates[[cov]])



      }

      inlabru::multiplot(plotlist = plotList)



    }
    else {


      ggplot() +
        boundaryComponent +
        meshComponent +
        speciesComponent



    }



  }
  ,
  addStructured = function(dataStructured, datasetType,
                            responseName, trialsName,
                            speciesName, coordinateNames) {

    if (missing(dataStructured)) stop('dataStructured needs to be provided')

    datasetName <- as.character(as.list(match.call())$dataStructured)

    message(paste('Adding dataset', datasetName, 'to the model.'))

    if (datasetName %in% names(private$dataStructured)) {

      warning('Dataset object already added to the model. Removing the previous version.')
      private$dataStructured[[datasetName]] <- NULL

      }

    if (!inherits(dataStructured,c('Spatial',
                  'data.frame','sf'))) stop('dataStructured needs to be either a SpatialPoints*, data.frame or sf object')

    if (missing(speciesName) && length(private$Species) > 1) stop('speciesName cannot be missing if the number of species in the model is greater than one.')

    if (missing(coordinateNames) && all(class(dataStructured) == 'data.frame')) stop('coordinateNames cannot be missing if dataStructured is a data.frame object.')
    else {

      if (inherits(dataStructured, 'Spatial')) coordinateNames <- colnames(dataStructured@coords)
      else coordinateNames <- colnames(st_coordinates(dataStructured))

    }

    if (is.null(private$Area)) stop('An area needs to be provided before adding species. This may be done with the `.$addArea` function.')

    if (missing(datasetType) || !datasetType %in% c('PO', 'PA', 'Counts')) stop('datasetType needs to be one of "PO", "PA" or "Counts".')

    if (missing(responseName) && datasetType %in% c('PA', 'Counts')) stop('responseName cannot be missing for PA and counts datasets.')

    if (datasetType == 'PA') responseNew <- private$responsePA
    if (datasetType == 'Counts') responseNew <- private$responseCounts
    if (datasetType == 'PO') responseNew <- NULL

    if (missing(trialsName)) trialsName <- NULL
    if (missing(responseName)) responseName <- NULL
    if (missing(speciesName)) speciesName <- NULL


    private$dataStructured[[datasetName]] <- formatStructured(data = dataStructured,
                                               type =  datasetType,
                                               varsOld = list(trials = trialsName,
                                                              response = responseName,
                                                              species = speciesName,
                                                              coordinates = coordinateNames),
                                               varsNew = list(coordinates = private$Coordinates,
                                                              response = responseNew,
                                                              trials = private$trialsName,
                                                              species = private$speciesName),
                                               projection = private$Projection,
                                               boundary = private$Area)



  }
  ,
  addMesh = function(Object,
                     cutoff,
                     max.edge,
                     offset,
                     ...) {

    if (missing(Object) &&
        is.null(Private$Area)) stop ('Please provide an inla.mesh object or use the other arguments to specift a mesh.')

    if (missing(Object)) {

      private$Mesh <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(as(private$Area, 'Spatial')),
                                         cutoff = cutoff,
                                         max.edge = max.edge,
                                         offset = offset,
                                         crs = private$Projection
                                          )

    }
    else {

      if (!inherits(Object, 'inla.mesh')) stop('Object provided is not an inla.mesh object.')

      private$Mesh <- Object

    }

    message('INLA mesh added successfully.')

    }
  ,
  addGBIF = function(Species = 'All',
                     gbifopts = list(), datasetName = NULL,
                     datasetType = 'PO',
                     responseCounts = 'individualCount', responsePA = 'occurrenceStatus',
                     assign2Global = FALSE, ...) {

  if (is.null(private$Area)) stop('An area needs to be provided before adding species. This may be done with the `.$addArea` function.')

  if (is.null(datasetName) && is.null(gbifopts$datasetKey)) stop('Please provide a name to give your dataset using datasetName. Alternatively, a name may be found using the `datasetKey` item in the gbifOptions argument.')

  if (Species == 'All') Species <- private$Species
  else if (!all(Species %in% private$Species)) stop ('Species provided not specified in startWorkflow().')

  if (datasetName %in% names(private$dataGBIF)) warning ('datasetName already provided before. The older dataset will therefore be removed.')

  ##To do here:
   #Keep only the relevant PA/Counts variables


  ##Change the speciesName argument to private$speciesName somewhere
   #Are we changing for Name or scientific name?
  private$dataGBIF[[datasetName]] <- obtainGBIF(query = Species,
                                                datasetName = datasetName,
                                                geometry = private$Area,
                                                gbifopts = gbifopts,
                                                projection = private$Projection,
                                                varsKeep = c(responseCounts, responsePA),
                                                ...)
  ##assign2Global and remove?
  if (assign2Global) assign(datasetName,  private$dataGBIF[[datasetName]], envir = globalenv())


  }
  ,
  addCovariates = function(Object = NULL,
                           worldClim = NULL,
                           Months = 'All',
                           Function = 'mean', ...) {

    if (is.null(private$Area)) stop("Area is required before obtaining covariates. Please use `.$addArea()`.")

    if (all(is.null(Object),
            is.null(worldClim))) stop ('One of object or worldClim is required.')


    if (length(worldClim) > 1) stop ('Please only add one worldClim variable at a time.')

    if (!worldClim %in% c("tavg", "tmin", "tmax",
                         "prec", "bio", "bioc",
                         "elev", "wind", "vapr", "srad")) stop('worldClim argument is not a valid option')


  if (is.null(Object)) {

    if (is.null(Months) && is.null(Function)) stop('Both of Months or Function need to be specified.')

    covDirectory <- paste0(private$Directory, '/Covariates')
    if (!file.exists(covDirectory)) dir.create(covDirectory)
    message(paste('Saved covariate objects may be found in', covDirectory))

    covRaster <- obtainCovariate(covariates = worldClim,
                                 countries = private$Countries,
                                 as.character(private$Projection),
                                 path = covDirectory)
    ##Crop this around boundary

      ##How does this work: Are a subset of months selected and then a function is applied?

      months <- c('January', 'February', 'March', 'April', 'May', 'June',
                  'July', 'August', 'September', 'October', 'November', 'December')

      if (!any(Months %in% c(months, 'All'))) stop ('Month provided is not valid.')

      if (all(Months == 'All')) covIndex <- 1:12
      else covIndex <- match(Months, months)

      private$Covariates[[worldClim]] <- do.call(Function, list(covRaster[[covIndex]]))


  }
    else {

      if (!inherits(Object, c('SpatRaster', 'SpatialPixelsDataFrame', 'raster'))) stop('Object needs to be either a SpatRaster, SpatialPixelsDataFrame or raster object.')

      Object <- as(Object, 'SpatRaster')

      Object <- terra::project(Object, private$Projection)

      #Check this for all classes
      maskedDF <- data.frame(mask(Object, private$Area))

      if (nrow(maskedDF) == 0) stop('Object provided does not mask the boundary.')
      #else get name of object and then save ti

    }

    }
  ,
  addArea = function(Object = NULL,
                     countryName = NULL) {

    if (all(is.null(Object),
            is.null(countryName))) stop ('One of object or countryName is required.')

    if (!is.null(private$Area)) {

      warning('Area already specified. Deleting all species occurance reccords added to the model.')

      private$dataStructured <- list()
      private$dataGBIF <- list()

    }

    if (!is.null(countryName)) private$Area <- obtainArea(name = countryName, projection =  private$Projection)
    else {

      if (!inherits(Object, 'Spatial') || !inherits(Object, 'sf')) stop('Object needs to be a sp or sf object.')

}

  }
  ,
  crossValidation = function(Method) {

    if (!all(Method %in% c('spatialBlock', 'Loo'))) stop('Output needs to be at least one of: spatialBlock, Loo.')

    private$CVMethod <- Method

  }
  ,
  modelOptions = function(ISDM = list(), INLA = list()) {

    if (!is.list(ISDM)) stop('ISDM needs to be a list of arguments to specify the model.')
    if (!is.list(INLA)) stop('INLA needs to be a list of INLA arguments to specify the model.')

    private$optionsISDM <- ISDM
    private$optionsINLA <- INLA

  }
  ,
  specifySpatial = function(sharedField, datasetName = list(), biasField = list()) {

    if (missing(sharedField) && length(datasetName) == 0 && length(biasField) == 0) stop('Please provide one of sharedField, datasetName or biasField.')

    if (!missing(sharedField)) {

      if (!inherits(sharedField, 'inla.spde2')) stop('sharedField needs to be a inla.spde2 object')
      private$sharedField <- sharedField


    }

    if (length(datasetName) > 0) {

      if (is.null(names(datasetName))) stop('datasetName is required to be a named list containing the names of the datasets to specify.')

      if (!all(names(datasetName) %in% c(names(private$dataGBIF), names(private$dataStructured)))) stop('datasetName contains a specification of a field for a dataset not present in the model.')



    }

    if (length(biasField) > 0) {

      if (is.null(names(datasetName))) stop('biasField is required to be a named list containing the names of the datasets to specify.')

      if (!all(names(biasField) %in% c(names(private$dataGBIF), names(private$dataStructured)))) stop('biasField contains a specification of a field for a dataset not present in the model.')

      }

  }
  ,
  biasFields = function(datasetName) {}
  ,
  workflowOutput = function(Output) {

  if (!all(Output %in% c('Model',
                         'Predictions',
                         'Maps',
                         'Cross-Validation'))) stop('Output needs to be at least one of: Model, Predictions, Maps or Cross-validation.')

    private$Output <- Output

  }
  ))

species_model$set('private', 'Area', NULL)
species_model$set('private', 'Covariates', NULL)
species_model$set('private', 'Mesh', NULL)
species_model$set('private', 'optionsISDM', list())
species_model$set('private', 'optionsINLA', list())
species_model$set('private', 'CVMethod', NULL)
species_model$set('private', 'Species', NULL)
species_model$set('private', 'Countries', NULL)
species_model$set('private', 'Output', NULL)
species_model$set('private', 'Projection', NULL)
species_model$set('private', 'dataStructured', list())
species_model$set('private', 'dataGBIF', list())
species_model$set('private', 'Quiet', TRUE)
species_model$set('private', 'Directory', getwd())
species_model$set('private', 'Project', NULL)
species_model$set('private', 'sharedField', NULL)
species_model$set('private', 'datasetName', list())


##Change all the variable names for the GBIF data too
species_model$set('private', 'Coordinates', c('Longitude', 'Latitude'))
species_model$set('private', 'responseCounts', 'individualCount')
species_model$set('private', 'responsePA', 'occurrenceStatus')
species_model$set('private', 'trialsName', 'numTrials')
species_model$set('private', 'speciesName', 'speciesName')

species_model$set('public', 'initialize', function(Countries, Species, nameProject,
                                                   Directory, Projection, Quiet = TRUE) {

  private$Projection <- Projection
  private$Quiet <- Quiet

  if (!missing(Countries)) {

    private$Countries <- Countries
    countryAdded <- self$addArea(countryName = Countries)

  }

  private$Species <- Species

  private$Directory <- Directory
  private$Project <- nameProject





})
