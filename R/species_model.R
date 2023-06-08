#' @title R6 class for creating a \code{species_model} object.
#' @description An object containing the data, covariates  and other relevant information to be used in the reproducible workflow. The function \link[intSDM]{startWorkflow} acts as a wrapper in creating one of these objects. This object has additional slot functions within, which allow for further specification and customization of the reproducible workflow.
#' @export
#' @importFrom R6 R6Class
#'
species_model <- R6::R6Class(classname = 'species_model', public = list(

#' @description initialize the species_model object.
#' @param Countries Name of the countries to include in the workflow.
#' @param Species Name of the species to include in the workflow.
#' @param nameProject Name of the project for the workflow.
#' @param Save Logical argument indicating if the model outputs should be saved.
#' @param Directory Directory where the model outputs should be saved.
#' @param Projection The coordinate reference system used in the workflow.
#' @param Quiet Logical variable indicating if the workflow should provide messages throughout the estimation procedure.
#'

  initialize = function(Countries, Species, nameProject, Save,
                        Directory, Projection, Quiet = TRUE) {

    private$Projection <- Projection
    private$Quiet <- Quiet

    if (!missing(Countries)) {

      #private$Countries <- Countries
      countryAdded <- self$addArea(countryName = Countries)

    }

    private$Species <- Species

    private$Directory <- Directory
    private$Project <- nameProject
    private$Save <- Save



  },

#' @description Prints the datasets, their data type and the number of observations, as well as the marks and their respective families.
#' @param ... Not used.
#' @import stats

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

      datasetNames <- unique(c(unlist(lapply(private$dataGBIF, function(x) names(x))),
                      unlist(lapply(private$dataStructured, function(x) names(x)))))

      cat(datasetNames, '\n\n') #DO an unlist



    }

    cat('Environmental covariates added: ')

    if (is.null(private$Covariates)) {

      cat('No covariates have been added to the workflow. Please add covariates using `.$addCovariates`.\n\n')

    } else {

      cat(names(private$Covariates), '\n\n')


    }

    cat('Model formula: ') ## Do this

    if (is.null(private$CVMethod)) cat('No Cross-validation specified. Please specify using `.$crossValidation`.\n\n')
    else cat('Cross-validation:', paste(private$CVMethod, collapse = ', '),'\n\n')

    if (is.null(private$Output)) cat('No Output has been specified. Please specify using `.$workflowOutput`.')
    else cat('Model output:', paste(private$Output, collapse = ', '))

   }
  ,

#' @description Makes a plot of the features used in the integrated model.
#' @param Mesh Add the mesh to the plot.
#' @param Boundary Add the boundary to the plot.
#' @param Species Add the species location data to the plot.
#' @param Covariates Add the spatial covariates to the plot.
#' @return A ggplot object.
#'
#' @import ggplot2
#' @import inlabru
#' @importFrom tidyterra geom_spatraster

  plot = function(Mesh = FALSE,
                  Boundary = TRUE,
                  Species = FALSE,
                  Covariates = FALSE) {

    if (sum(Mesh, Boundary, Species, Covariates) == 0) stop('At least one of Mesh, Boundary, Species or Covariates needs to be chosen.')

    if (Mesh) {

      if (is.null(private$Mesh)) stop('No mesh object provided. Please add one using `.$addMesh()`.')

      meshComponent <- inlabru::gg(private$Mesh)

    } else meshComponent <- NULL

    if (Boundary) {

      if (is.null(private$Area)) stop('No boundary object proived. Please add one using `.$addArea()`.')

      boundaryComponent <- geom_sf(data = sf::st_boundary(private$Area))

    } else boundaryComponent <- NULL

    if (Species) {

      if (all(is.null(private$dataGBIF), is.null(private$dataStructured))) stop('No data objects provided. Please add using either `.$addGBIF` or `.$addStructured`.')

      spatData <- vector(mode = 'list', length = length(private$Species))
      names(spatData) <-  sub(' ', '_', private$Species)

      for (species in sub(' ', '_', private$Species)) {

        speData <- lapply(append(private$dataGBIF[[species]],
                                 private$dataStructured[[species]]), function(x) x$geometry)

        namesspeData <- names(speData)
        namesTimes <- rep(namesspeData, times = unlist(lapply(speData, length)))

        spatData[[species]] <- sf::st_as_sf(do.call(c, speData))

        spatData[[species]]$.__species_index_var <- species
        spatData[[species]]$.__names_index_var <- namesTimes

      }

      plotSpecies <- do.call(rbind, spatData)
      speciesComponent <- geom_sf(data = plotSpecies, aes(col = .__names_index_var))
      guidesComponent <- guides(col = guide_legend(title="Dataset Name")) ##Need to convert this speciesName arg for all species
      facetComponent <- facet_wrap( ~ .__species_index_var)
      }

    else {

        speciesComponent <- NULL
        guidesComponent <- NULL
        facetComponent <- NULL

      }

    if (Covariates) {

      if (is.null(private$Covariates)) stop('No covariates object provided. Please add using `.$addCovariates()`.')

      plotList <- list()

      for (cov in names(private$Covariates)) {


        plotList[[cov]] <- ggplot() +
                           tidyterra::geom_spatraster(data = private$Covariates[[cov]]) +
                           boundaryComponent +
                           meshComponent +
                           speciesComponent +
                           guidesComponent



      }

      if (length(plotList) > 1) inlabru::multiplot(plotlist = plotList)
      else return(plotList[[1]])



    }
    else {


      plotStr <- ggplot() +
                 boundaryComponent +
                 meshComponent +
                 speciesComponent +
                 guidesComponent +
                 facetComponent

      return(plotStr)


      }



  }
  ,

#' @description The function is used to convert structured datasets into a framework which is usable by the model. The three types of structured data allowed by this function are present absence (PA), present only (PO) and counts/abundance datasets, which are controlled using the \code{datasetType} argument. The other arguments of this function are used to specify the appropriate variable (such as response name, trial name, species name and coordinate name) names in these datasets.
#'
#' @param dataStructured The dataset used in the model. May be either a \code{data.frame}, \code{sf} or \code{SpatialPoints*} object.
#' @param datasetType A vector which gives the type of dataset. Must be either \code{'count'}, \code{'PO'} or \code{'PA'}.
#' @param responseName Name of the response variable in the dataset. If \code{dataType} is \code{'PO'}, then this argument may be missing.
#' @param trialsName Name of the trial name variable in the \code{PA} datasets.
#' @param speciesName Name of the species variable name in the datasets.
#' @param coordinateNames Names of the coordinate vector in the dataset.
#'
#' @import methods
#' @import sp
#' @import sf
#'
  addStructured = function(dataStructured, datasetType,
                            responseName, trialsName,
                            speciesName, coordinateNames) {

    if (missing(dataStructured)) stop('dataStructured needs to be provided')

    datasetName <- as.character(as.list(match.call())$dataStructured)

    if (!private$Quiet) message(paste('Adding dataset', datasetName, 'to the model.'))

    if (datasetName %in% names(private$dataStructured)) {

      warning('Dataset object already added to the model. Removing the previous version.')
      private$dataStructured[[datasetName]] <- NULL

      }

    if (!inherits(dataStructured,c('Spatial',
                  'data.frame','sf'))) stop('dataStructured needs to be either a SpatialPoints*, data.frame or sf object')

    if (missing(speciesName)) stop('speciesName cannot be missing.')

    if (missing(coordinateNames) && all(class(dataStructured) == 'data.frame')) stop('coordinateNames cannot be missing if dataStructured is a data.frame object.')
    else {

      if (inherits(dataStructured, 'Spatial')) coordinateNames <- colnames(dataStructured@coords)
      if (inherits(dataStructured, 'sf')) coordinateNames <- colnames(st_coordinates(dataStructured))

    }

    if (is.null(private$Area)) stop('An area needs to be provided before adding species. This may be done with the `.$addArea` function.')

    if (missing(datasetType) || !datasetType %in% c('PO', 'PA', 'Counts')) stop('datasetType needs to be one of "PO", "PA" or "Counts".')

    if (missing(responseName) && datasetType %in% c('PA', 'Counts')) stop('responseName cannot be missing for PA and counts datasets.')

    if (datasetType == 'PA') responseNew <- private$responsePA
    if (datasetType == 'Counts') responseNew <- private$responseCounts
    if (datasetType == 'PO') responseNew <- '.__PORESP.__' #NULL

    if (missing(trialsName)) trialsName <- NULL
    if (missing(responseName)) responseName <- NULL
    #if (missing(speciesName)) speciesName <- NULL

    dataStructured[[speciesName]] <- sub(" ", "_", dataStructured[[speciesName]])

    uniqueSpecies <- unique(dataStructured[[speciesName]])

    if (!all(uniqueSpecies %in% sub(" ", "_", private$Species))) {

      warning('Species found in dataset not specified in the original startWorkflow call. Removing observations for those species.')

      dataStructured <- dataStructured[dataStructured[[speciesName]] %in% sub(" ", "_", private$Species), ]

      uniqueSpecies <- unique(dataStructured[[speciesName]])

      if (nrow(dataStructured) == 0) stop('All species removed from this dataset. Please specify all the species name using the "Species" argument in startWorkflow.')

    }


  for (species in uniqueSpecies) {


    private$dataStructured[[species]][[datasetName]] <- formatStructured(data = dataStructured,
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

  private$datasetName <- c(datasetName, private$datasetName)


  }
  ,

#' @description Function to add an \code{inla.mesh} object to the workflow. The user may either add their own mesh to the workflow, or use the arguments of this function to help create one.
#' @param Object An \code{inla.mesh} object to add to the workflow.
#' @param ... Additional arguments to pass to \code{INLA}'s \code{inla.mesh.2d}. Use \code{?inla.mesh.2d} to find out more about the different arguments.
#'
#' @import INLA
#'

  addMesh = function(Object,
                     ...) {

    if (missing(Object) &&
        is.null(private$Area)) stop ('Please provide an inla.mesh object or use the ... argument to specify a mesh.')

    if (missing(Object)) {

      meshArgs <- list(...)
      if (length(meshArgs) == 0) stop('Please provide ... to specify the mesh construction. See ?inla.mesh.2d for more details.')

      meshObj <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(as(private$Area[1], 'Spatial')),
                                    crs = private$Projection,
                                         ...
                                          )

      private$Mesh <- meshObj

    }
    else {

      if (!inherits(Object, 'inla.mesh')) stop('Object provided is not an inla.mesh object.')

      private$Mesh <- Object

    }

    if (!private$Quiet) message('INLA mesh added successfully.')

    }
  ,

#' @description Function to add species occurrence records from GBIF (using the \code{rgbif} package) to the reproducible workflow. The arguments for this function are used to either filter the GBIF records, or to specify the characteristics of the observation model.
#' @param Species The names of the species to include in the workflow (initially specified using \link[intSDM]{startWorkflow}). Defaults to \code{All}, which will find occurrence records for all specie specified in \link[intSDM]{startWorkflow}.
#' @param datasetName The name to give the dataset obtained from GBIF. Cannot be \code{NULL}.
#' @param datasetType The data type of the dataset. Defaults to \code{PO}, but may also be \code{PA} or \code{Counts}.
#' @param responseCounts Name of the response variable for the counts data. Defaults to the standard Darwin core value \code{individualCounts}.
#' @param responsePA Name of the response variable for the PA data. Defaults to the standard Darwin core value \code{occurrenceStatus}.
#' @param assign2Global Assign the dataset to the global environment. The object will be assigned to an object specified using the \code{datasetName} object.
#' @param ... Additional arguments to specify the \link[rgbif]{occ_data} function from \code{rgbif}. See \code{?occ_data} for more details.

addGBIF = function(Species = 'All', datasetName = NULL,
                   datasetType = 'PO',
                   responseCounts = 'individualCount', responsePA = 'occurrenceStatus',
                   assign2Global = FALSE, ...) {

  if (is.null(private$Area)) stop('An area needs to be provided before adding species. This may be done with the `.$addArea` function.')

  if (is.null(datasetName)) stop('Please provide a name to give your dataset using datasetName.')

  if (Species == 'All') Species <- private$Species
  else if (!all(Species %in% private$Species)) stop ('Species provided not specified in startWorkflow().')

  if (datasetName %in% names(private$dataGBIF)) warning ('datasetName already provided before. The older dataset will therefore be removed.')

  ##To do here:
   #Keep only the relevant PA/Counts variables

  if (!datasetType %in% c('PO', 'PA', 'Counts')) stop('datasetType needs to be one of PO, PA and Counts.')

  if (datasetType == 'PA' && is.null(responsePA)) stop('responsePA cannot be NULL if datasetType is PA.')

  if (datasetType == 'Counts' && is.null(responseCounts)) stop('responseCounts cannot be NULL if datasetType is Counts.')


  ##Change the speciesName argument to private$speciesName somewhere
   #Are we changing for Name or scientific name?

  for (speciesName in Species) {

  if (!private$Quiet) message(paste('Finding GBIF observations for:', speciesName,'\n'))

  GBIFspecies <- obtainGBIF(query = speciesName,
                                                #datasetName = datasetName,
                            geometry = private$Area,
                                                #country = private$Countries,
                            projection = private$Projection,
                                                #varsKeep = c(responseCounts, responsePA),
                            datasettype = datasetType,
                            ...)

  if (length(private$dataGBIF[[sub(" ", '_', speciesName)]]) > 0) {

    anySame <- st_equals_exact(do.call(c, lapply(private$dataGBIF[[sub(" ", '_', speciesName)]], function(x) st_geometry(x))),
                               GBIFspecies,
                               par = 1e-3)

    if (!identical(unlist(anySame), 'integer(0)')) {

      warning('Removing duplicate observations obtained from previous calls of `.$addGBIF`')

      GBIFspecies <- GBIFspecies[-c(unlist(anySame)),]


    }


  }

  if (!assign2Global) {

  if (nrow(GBIFspecies) == 0) warning('All species observations were removed due to duplicates')
  else {

    private$dataGBIF[[sub(" ", '_', speciesName)]][[datasetName]] <- GBIFspecies
    private$classGBIF[[sub(" ", '_', speciesName)]][[datasetName]] <- datasetType

  }
  } else assign(datasetName,  private$dataGBIF[[speciesName]][[datasetName]], envir = globalenv())

  }

  private$datasetName <- c(datasetName, private$datasetName)

  }
  ,

#' @description Function to add spatial covariates to the workflow. The covariates may either be specified by the user, or they may come from worldClim obtained with the \code{geodata} package.
#' @param Object A object of class: \code{spatRaster}, \code{SpatialPixelsDataFrame} or \code{raster} containing covariate information across the area. Note that this function will check if the covariates span the boundary area, so it may be preferable to add your own boundary using \code{`.$addArea`} if this argument is specified.
#' @param worldClim Name of the worldClim to include in the model. See \code{?worldclim_country} from the \code{geodata} package for more information.
#' @param Months The months to include the covariate for. Defaults to \code{All} which includes covariate layers for all months.
#' @param Function The function to aggregate the temporal data into one layer. Defaults to \code{mean}.
#' @param ... Not used.
#'
#' @import terra
#'
  addCovariates = function(Object = NULL,
                           worldClim = NULL,
                           Months = 'All',
                           Function = 'mean', ...) {

    if (is.null(private$Area)) stop("Area is required before obtaining covariates. Please use `.$addArea()`.")

    if (all(is.null(Object),
            is.null(worldClim))) stop ('One of object or worldClim is required..')


    if (length(worldClim) > 1) stop ('Please only add one worldClim variable at a time.')


  if (is.null(Object)) {

    if (!worldClim %in% c("tavg", "tmin", "tmax",
                          "prec", "bio", "bioc",
                          "elev", "wind", "vapr", "srad")) stop('worldClim argument is not a valid option.')

    if (is.null(Months) && is.null(Function)) stop('Both of Months or Function need to be specified.')

    months <- c('January', 'February', 'March', 'April', 'May', 'June',
                'July', 'August', 'September', 'October', 'November', 'December')

    if (!any(Months %in% c(months, 'All'))) stop ('Month provided is not valid.')

    if (all(Months == 'All')) covIndex <- 1:12
    else covIndex <- match(Months, months)

    covDirectory <- paste0(private$Directory, '/Covariates')
    dir.create(covDirectory, showWarnings = FALSE)
    if (!dir.exists(covDirectory)) covDirectory <- getwd()#dir.create(covDirectory)
    if (!private$Quiet) message(paste('Saved covariate objects may be found in', covDirectory))

    if (is.null(private$Countries)) stop('Please specify a country first before obtaining a covariate layer. This may be done using either startWorkflow or through `.$addArea`.')

    covRaster <- obtainCovariate(covariates = worldClim,
                                 countries = private$Countries,
                                 as.character(private$Projection),
                                 path = covDirectory)

    covRaster <- terra::mask(covRaster, private$Area)

    covRaster <- terra::app(covRaster[[covIndex]], fun = Function)
    names(covRaster) <- worldClim

    private$Covariates[[worldClim]] <- covRaster


  }
    else {

      if (!inherits(Object, c('SpatRaster', 'Spatial', 'Raster'))) stop('Object needs to be either a SpatRaster, SpatialPixelsDataFrame or raster object.')

      if (!inherits(Object, 'SpatRaster')) Object <- as(Object, 'SpatRaster')

      Object <- terra::project(Object, private$Projection)

      if (length(names(Object)) > 1) stop('Please provide each covariate into the workflow as their own object.')

      #Check this for all classes
      maskedDF <- terra::mask(Object, private$Area)

      if (all(is.na(terra::values(maskedDF)))) stop('The covariate provided and the area specified do not match.')

      private$Covariates[[names(Object)]] <- Object
      #else get name of object and then save ti

    }

    }
  ,

#' @description Function to add a boundary around the study area. This function allows the user to either add their own boundary object, or obtain a country's boundary using \link[giscoR]{gisco_get_countries} from the \code{giscoR} package.
#' @param Object A \code{sf} or \code{SpatialPolygons} object of the boundary surrounding the study area.
#' @param countryName Name of the countries to obtain a boundary for. This argument will then use the \link[giscoR]{gisco_get_countries} function from the \code{giscoR} package to obtain a boundary.
#' @param ... Additional arguments passed to \link[giscoR]{gisco_get_countries}.
#'
#' @import sf
#'
  addArea = function(Object = NULL,
                     countryName = NULL,
                     ...) {

    if (all(is.null(Object),
            is.null(countryName))) stop ('One of object or countryName is required.')

    if (!is.null(private$Area)) {

      warning('Area already specified. Deleting all species occurance reccords added to the model.')

      private$dataStructured <- list()
      private$dataGBIF <- list()

    }

    if (!is.null(countryName)) {

      private$Area <- obtainArea(name = countryName, projection =  private$Projection, ...)
      private$Countries <- countryName

      }
    else {

      if (!inherits(Object, 'Spatial') && !inherits(Object, 'sf')) stop('Object needs to be a sp or sf object.')
      if (inherits(Object, 'Spatial')) Object <- as(Object, 'sf')

      Object <- sf::st_transform(Object, as.character(private$Projection))

      private$Area <- Object

    }

    if (!private$Quiet) message('Boundry object added successfully.')

  }
  ,

#' @description Function to add a spatial cross validation method to the workflow.
#' @param Method The spatial cross-validation methods to use in the workflow. May be at least one of \code{spatialBlock} or \code{Loo} (leave-one-out). See the \code{PointedSDMs} package for more details.
#' @param blockOptions A list of options to specify the spatial block cross-validation. Must be a named list with arguments specified for: \code{k}, \code{rows_cols}, \code{plot}, \code{seed}. See \code{blockCV::cv_spatial} for more information.
#'
#' @import blockCV
#'
  crossValidation = function(Method, blockOptions = list(k = 5, rows_cols = c(4,4), plot = FALSE, seed = NULL)) {

    if (!all(Method %in% c('spatialBlock', 'Loo'))) stop('Output needs to be at least one of: spatialBlock, Loo.')

    if ('spatialBlock' %in% Method) {


     if (is.null(blockOptions$seed)) blockOptions$seed <- round(abs(rnorm(n = 1, mean = 100000, sd = 100000)))

     if (is.null(blockOptions$k) || is.null(blockOptions$rows_cols)) stop('Please provide both k and rows_cols in blockOptions.')
     if (is.null(blockOptions$plot)) blockOptions$plot <- FALSE

     private$blockOptions <- blockOptions

     if (blockOptions$plot) {

     boundData <- vector(mode = 'list', length =  length(private$Species))
     spatPlot <- vector(mode = 'list', length =  length(private$Species))
     names(spatPlot) <- paste('Printing plot for',private$Species)


       for (species in 1:length(private$Species)) {

       spatData <- lapply(append(private$dataGBIF[[species]], private$dataStructured[[species]]), function(x) x$geometry)

       boundData <- st_as_sf(do.call(c, spatData))


       blocks <- R.devices::suppressGraphics(blockCV::cv_spatial(x = boundData,
                                                                 k = blockOptions$k, rows_cols = blockOptions$rows_cols,
                                                                 progress = FALSE, seed = blockOptions$seed,
                                                                 report = FALSE, plot = TRUE))

       boundData$.__block_index <- blocks$folds_ids

       blocksPlot <- blockCV::cv_plot(blocks)

        spatPlot[[species]] <-  ggplot() +
         geom_sf(data = blocksPlot$data, aes(col = as.character(folds)), size = 2) +
         geom_sf_text(data = blocksPlot$data, aes(label = folds)) +
         geom_sf(data = sf::st_boundary(private$Area)) +
         geom_sf(data = boundData, aes(col = as.character(.__block_index))) +
         ggtitle(paste('Cross-validation blocking for', private$Species[species])) +
         guides(col=guide_legend(title="Folds"))



       }


     return(spatPlot)


     }


    }

    private$CVMethod <- Method

  }
  ,

#' @description Function to specify model options for the \code{INLA} and \code{PointedSDMs} parts of the model.
#' @param ISDM Arguments to specify in \link[PointedSDMs]{intModel} from the \code{PointedSDMs} function. This argument needs to be a named list of the following options: \code{pointCovariates}, \code{pointsIntercept}, \code{pointsSpatial} or \code{copyModel}. See \code{?intModel} for more details.
#' @param INLA Options to specify in \link[INLA]{inla} from the \code{INLA} function. See \code{?inla} for more details.
#'
  modelOptions = function(ISDM = list(),
                          INLA = list()) {

    if (!is.list(ISDM)) stop('ISDM needs to be a list of arguments to specify the model.')

    if (any(!names(ISDM) %in% c('pointCovariates', 'pointsIntercept', #Remove pointCovariates perhaps?
                                'pointsSpatial', 'copyModel'))) stop('ISDM needs to be a named list with at least one of the following options: "pointCovariates", "pointsIntercept", "pointsSpatial" or "copyModel".')

    if (!is.list(INLA)) stop('INLA needs to be a list of INLA arguments to specify the model.')

    if (length(ISDM) > 0) private$optionsISDM <- ISDM
    if (length(INLA) > 0) private$optionsINLA <- INLA

  }
  ,

#' @description Function to specify pc priors for the shared random field in the model. See \code{?INLA::inla.spde2.pcmatern} for more details.
#' @param ... Arguments passed on to \link[INLA]{inla.spde2.pcmatern}.
#'
#' @import INLA
  specifySpatial = function(...) {


    anyIn <- list(...)
    if (length(anyIn) == 0) stop('Please provide arguments to customize the INLA spde object using the ... argument.')

    if (is.null(private$Mesh)) stop('Please add an INLA mesh before customizing the spatial fields. This may be done with the `.$addMesh` function.')

    private$sharedField <- INLA::inla.spde2.pcmatern(mesh = private$Mesh, ...)

  }
  ,

#' @description Function to add bias fields to the model.
#' @param datasetName Name of the dataset to add a bias field to.
#' @param ... Additional arguments passed on to \link[INLA]{inla.spde2.pcmatern} to customize the priors for the pc matern for the bias fields.
#'
#' @import INLA

  biasFields = function(datasetName, ...) {

    if (!all(datasetName %in% private$datasetName)) stop('Dataset specified for bias field not included in the workflow.')

    private$biasNames <- unique(c(datasetName, private$biasNames))

    if (length(list(...)) > 0) {

      if (is.null(private$Mesh)) stop('Please add an INLA mesh before customizing the spatial fields. This may be done with the `.$addMesh` function.')

      biasModels <- INLA::inla.spde2.pcmatern(mesh = private$Mesh, ...)

      for (dataset in datasetName) {

      private$biasFieldsSpecify[[dataset]] <- biasModels


      }

    }


  }
  ,

#' @description Function to specify the workflow output from the model. This argument must be at least one of: \code{'Model'}, \code{'Prediction'}, \code{'Maps'} and \code{'Cross-validation'}.
#' @param Output The names of the outputs to give in the workflow. Must be at least one of: \code{'Model'}, \code{'Prediction'}, \code{'Maps'} and \code{'Cross-validation'}.
#'
  workflowOutput = function(Output) {

  if (!all(Output %in% c('Model',
                         'Predictions',
                         'Maps',
                         'Cross-validation'))) stop('Output needs to be at least one of: Model, Predictions, Maps or Cross-validation.')

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
#species_model$set('private', 'datasetFieldsSpecify', list())
species_model$set('private', 'biasFieldsSpecify', list())
species_model$set('private', 'datasetName', NULL)
species_model$set('private', 'Save', TRUE)
species_model$set('private', 'biasNames', NULL)
species_model$set('private', 'blockOptions', list())
species_model$set('private', 'classGBIF', list())



##Change all the variable names for the GBIF data too
species_model$set('private', 'Coordinates', c('Longitude', 'Latitude'))
species_model$set('private', 'responseCounts', 'individualCount')
species_model$set('private', 'responsePA', 'occurrenceStatus')
species_model$set('private', 'trialsName', 'numTrials')
species_model$set('private', 'speciesName', 'speciesName')

