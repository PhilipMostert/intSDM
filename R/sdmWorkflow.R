#' @title \code{sdmWorkflow}: Function to compile the reproducible workflow.
#' @description This function is used to compile the reproducible workflow from the \code{R6} object created with \code{startFunction}. Depending on what was specified before, this function will estimate the integrated species distribution model, perform cross-validation, create predictions from the model and plot these predictions.
#' @param Workflow The \code{R6} object created from \code{startWorkflow}. This object should contain all the data and model information required to estimate and specify the model.
#'
#' @import PointedSDMs
#'
#' @return The return of the function depends on the argument \code{Save} from the \code{startWorkflow} function. If this argument is \code{FALSE} then the objects will be saved to the specidfied directory. If this argument is \code{TRUE} then a list of different outcomes from the workflow will be returned.
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace('INLA')) {
#'
#' workflow <- startWorkflow(Species = 'Fraxinus excelsior',
#'                           Projection = '+proj=longlat +ellps=WGS84',
#'                           Save = FALSE,
#'                           saveOptions = list(projectName = 'example'))
#' workflow$addArea(countryName = 'Sweden')
#'
#' workflow$addGBIF(datasetName = 'exampleGBIF',
#'                  datasetType = 'PA',
#'                  limit = 10000,
#'                  coordinateUncertaintyInMeters = '0,50')
#' workflow$addMesh(cutoff = 20000,
#'                  max.edge=c(60000, 80000),
#'                  offset= 100000)
#' workflow$workflowOutput('Model')
#'
#' Model <- sdmWorkflow(workflow)
#'
#' }
#' }

sdmWorkflow <- function(Workflow = NULL) {

  modDirectory <- Workflow$.__enclos_env__$private$Directory
  saveObjects <- Workflow$.__enclos_env__$private$Save
  Quiet <- Workflow$.__enclos_env__$private$Quiet

  if (is.null(Workflow$.__enclos_env__$private$Mesh)) stop('An inla.mesh object is required before any analysis is completed. Please add using the `.$addMesh` function.')
  if (is.null(Workflow$.__enclos_env__$private$Output)) stop('A model output needs to be specified before any analysis is completed. Please add using the `.$workflowOutput` function.')

  Oputs <- Workflow$.__enclos_env__$private$Output
  outputList <-list()

  if (is.null(Workflow$.__enclos_env__$private$CVMethod) && 'Cross-Validation' %in% Workflow$.__enclos_env__$private$Output) stop('Cross-validation specified as model output but no method provided. Please specify cross-validation method using the `.$crossValidation` function.')

  ##Need to do it on a species by species basis:

  if (length(Workflow$.__enclos_env__$private$optionsISDM) > 0) {

    if (!is.null(Workflow$.__enclos_env__$private$optionsISDM[['pointsSpatial']])) .__pointsSpatial.__ <- Workflow$.__enclos_env__$private$optionsISDM$pointsSpatial
    else .__pointsSpatial.__ <- 'copy'#'shared' #Should this be copy?

    if (!is.null(Workflow$.__enclos_env__$private$optionsISDM[['copyModel']])) .__copyModel.__ <- eval(parse(text = Workflow$.__enclos_env__$private$optionsISDM['copyModel']))
    else .__copyModel.__ <- eval(parse(text = 'list(beta = list(fixed = FALSE))'))

    if (!is.null(Workflow$.__enclos_env$private$optionsISDM[['pointsIntercept']])) .__pointsIntercept.__ <- Workflow$.__enclos_env$private$optionsISDM$pointsIntercept
    else .__pointsIntercept.__ <- TRUE


  } else {

    .__pointsSpatial.__ <- 'copy'
    .__copyModel.__ <- eval(parse(text = 'list(beta = list(fixed = FALSE))'))
    .__pointsIntercept.__ <- TRUE

  }

  .__mesh.__ <- Workflow$.__enclos_env__$private$Mesh
  .__proj.__ <- Workflow$.__enclos_env__$private$Projection
  .__coordinates.__ <- Workflow$.__enclos_env__$private$Coordinates
  .__responsePA.__ <- Workflow$.__enclos_env__$private$responsePA
  .__responseCounts.__ <- Workflow$.__enclos_env__$private$responseCounts
  .__trialsName.__ <- Workflow$.__enclos_env__$private$trialsName

  if (!is.null(Workflow$.__enclos_env__$private$Covariates)) spatCovs <- terra::rast(Workflow$.__enclos_env__$private$Covariates)
  else spatCovs <- NULL

  for (species in unique(c(names(Workflow$.__enclos_env__$private$dataGBIF),
                           names(Workflow$.__enclos_env__$private$dataStructured)))) {

   speciesNameInd <- sub(' ', '_', species)
   if (saveObjects) dir.create(path = paste0(modDirectory, '/', speciesNameInd))

   if (!Quiet) message(paste('\nStarting model for', species, '\n\n'))

   speciesDataset <- append(Workflow$.__enclos_env__$private$dataGBIF[[species]],
                            Workflow$.__enclos_env__$private$dataStructured[[species]])



  if (length(speciesDataset) == 0) stop('No data added to the model. Please add data using `.$addGBIF` or `.$addStructured`.')

   if (!Quiet) message('\nInitializing model', '\n\n')

   if (length(speciesDataset) == 1) {

     initializeModel <- dataSDM$new(coordinates = .__coordinates.__, projection = .__proj.__,
                            Inlamesh = .__mesh.__, initialnames = names(speciesDataset),
                            responsecounts = .__responseCounts.__,
                            responsepa = .__responsePA.__,
                            marksnames = NULL,
                            marksfamily = NULL,
                            pointcovariates = NULL,
                            trialspa = .__trialsName.__,
                            trialsmarks = NULL,
                            spatial = .__pointsSpatial.__,
                            marksspatial = NULL,
                            speciesname = NULL,
                            intercepts = .__pointsIntercept.__,
                            marksintercepts = NULL,
                            spatialcovariates = spatCovs,
                            boundary = NULL,
                            ips = NULL,
                            temporal = NULL,
                            temporalmodel = NULL,
                            speciesspatial = NULL,
                            offset = NULL,
                            copymodel = .__copyModel.__)
     initializeModel$addData(speciesDataset)
   }
   else {

  initializeModel <- PointedSDMs::intModel(speciesDataset, Mesh = .__mesh.__, Projection = .__proj.__, Coordinates = .__coordinates.__,
                                            responsePA = .__responsePA.__, responseCounts = .__responseCounts.__,
                                            trialsPA = .__trialsName.__, pointsSpatial = .__pointsSpatial.__,
                                            pointsIntercept = .__pointsIntercept.__ ,
                                            copyModel = .__copyModel.__,
                                            spatialCovariates = spatCovs)

   }

  if (!is.null(Workflow$.__enclos_env__$private$sharedField)) initializeModel$spatialFields$sharedField$sharedField <- Workflow$.__enclos_env__$private$sharedField

  if (!is.null(Workflow$.__enclos_env__$private$biasNames)) {

    initializeModel$addBias(datasetNames = Workflow$.__enclos_env__$private$biasNames)

    if (!is.null(Workflow$.__enclos_env__$private$biasFieldsSpecify)) {

      for (biasName in names(Workflow$.__enclos_env__$private$biasFieldsSpecify)) {

        initializeModel$spatialFields$biasFields[[biasName]] <- Workflow$.__enclos_env__$private$biasFieldsSpecify[[biasName]]

        }


    }



  }

  if ('Cross-validation' %in% Oputs && 'spatialBlock' %in%Workflow$.__enclos_env__$private$CVMethod) {

    initializeModel$spatialBlock(k = Workflow$.__enclos_env__$private$blockOptions$k,
                                 rows_cols = Workflow$.__enclos_env__$private$blockOptions$rows_cols,
                                 seed = Workflow$.__enclos_env__$private$blockOptions$seed)

  }

  message('\nEstimating ISDM:\n\n')

  PSDMsMOdel <- PointedSDMs::fitISDM(initializeModel,
                                     options = Workflow$.__enclos_env__$private$optionsINLA)

  if ('Model' %in% Oputs) {

    if (saveObjects) {

    if (!Quiet) message('\nSaving Model object:', '\n\n')
    saveRDS(object = PSDMsMOdel, file = paste0(modDirectory,'/', speciesNameInd, '/intModel.rds'))

    } else outputList[[speciesNameInd]][['Model']] <- PSDMsMOdel

  }

  if ('Richness' %in% Oputs) {

    message('Richness output not available yet!')

    #message('\nCreating richness model:\n\n')

    #message('\nCreating richness maps:\n\n')

    #richnessMap <- predict(richnessModel)

  }

  if ('Cross-validation' %in% Oputs) {


    if ('spatialBlock' %in% Workflow$.__enclos_env__$private$CVMethod) {

      if (!Quiet) message('\nEstimating spatial block cross-validation:\n\n')
      spatialBlockCV <- PointedSDMs::blockedCV(initializeModel)

      if (saveObjects) {

      if (!Quiet) message('\nSaving spatial blocked cross-validation object:', '\n\n')
      saveRDS(object = spatialBlockCV, file = paste0(modDirectory,'/', speciesNameInd, '/spatialBlock.rds'))

      } else outputList[[speciesNameInd]][['spatialBlock']] <- spatialBlockCV

      rm(spatialBlockCV)
    }

    if ('Loo' %in% Workflow$.__enclos_env__$private$CVMethod) {

      if (!Quiet) message('\nEstimating leave-one-out cross-validation:\n\n')
      LooCV <- PointedSDMs::datasetOut(model = PSDMsMOdel)

      if (saveObjects) {

      if (!Quiet) message('\nSaving leave-one-out cross-validation object:', '\n\n')
      saveRDS(object = LooCV, file = paste0(modDirectory,'/', speciesNameInd, '/LooCV.rds'))

      } else outputList[[speciesNameInd]][['LooCV']] <- LooCV

      rm(LooCV)

    }

  }

    if (any(c('Predictions', 'Maps') %in% Oputs)) {

      if (!Quiet) message('\nPredicting model:\n\n')
      .__mask.__ <- as(Workflow$.__enclos_env__$private$Area, 'Spatial')
      Predictions <- predict(PSDMsMOdel, data = inlabru::fm_pixels(mesh = .__mesh.__,
                                                                mask = .__mask.__),
                             predictor = TRUE)

      if (saveObjects) {

      if (!Quiet)  message('\nSaving predictions object:', '\n\n')
      saveRDS(object = Predictions, file = paste0(modDirectory,'/', speciesNameInd, '/Predictions.rds'))

      } else outputList[[speciesNameInd]][['Predictions']] <- Predictions

    }

    if ('Maps' %in% Oputs) {

      if (!Quiet) message("\nPlotting Maps:\n\n")

      if (saveObjects) {

      if (!Quiet) message('\nSaving plots object:', '\n\n')
       plot(Predictions)
       ggsave(filename = paste0(modDirectory,'/', speciesNameInd, '/Map.png'))

      }
      else outputList[[speciesNameInd]][['Maps']] <- plot(Predictions, plot = FALSE)

      rm(Predictions)


    }

  }

  if (!saveObjects) return(outputList)


  }
