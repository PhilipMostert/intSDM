#' @title \code{sdmWorkflow}: Function to compile the reproducible workflow.
#' @description This function is used to compile the reproducible workflow from the \code{R6} object created with \code{startFunction}. Depending on what was specified before, this function will estimate the integrated species distribution model, perform cross-validation, create predictions from the model and plot these predictions.
#' @param Workflow The \code{R6} object created from \code{startWorkflow}. This object should contain all the data and model information required to estimate and specify the model.
#' @param predictionDim The pixel dimensions for the prediction maps. Defaults to \code{c(150, 150)}.
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

sdmWorkflow <- function(Workflow = NULL,
                        predictionDim = c(150, 150)) {

  modDirectory <- Workflow$.__enclos_env__$private$Directory
  saveObjects <- Workflow$.__enclos_env__$private$Save
  Quiet <- Workflow$.__enclos_env__$private$Quiet

  if (is.null(Workflow$.__enclos_env__$private$Mesh)) stop('An inla.mesh object is required before any analysis is completed. Please add using the `.$addMesh` function.')
  if (is.null(Workflow$.__enclos_env__$private$Output)) stop('A model output needs to be specified before any analysis is completed. Please add using the `.$workflowOutput` function.')

  Oputs <- Workflow$.__enclos_env__$private$Output
  outputList <-list()

  if ('Bias' %in% Oputs && is.null(Workflow$.__enclos_env__$private$biasNames)) stop('Bias specified as output but no bias fields were added. Please do this with .$biasFields')

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
                            speciesintercept = FALSE,
                            speciesindependent = FALSE,
                            speciesenvironment = FALSE,
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

    if (!any(names(speciesDataset) %in% Workflow$.__enclos_env__$private$biasNames)) {

      warning('Bias fields specified for datasets not in the model. Turning bias off.')

      biasIn <- FALSE

      }

    else {

    biasIn <- TRUE

    biasSubset <- Workflow$.__enclos_env__$private$biasNames[Workflow$.__enclos_env__$private$biasNames %in% names(speciesDataset)]

    initializeModel$addBias(datasetNames = biasSubset, copyModel = Workflow$.__enclos_env__$private$biasFieldsCopy, shareModel = Workflow$.__enclos_env__$private$biasFieldsShare)

    if (!is.null(Workflow$.__enclos_env__$private$biasFieldsSpecify)) {

      if (any(names(speciesDataset) %in% names(Workflow$.__enclos_env__$private$biasFieldsSpecify))) {

        specifySubset <- names(Workflow$.__enclos_env__$private$biasFieldsSpecify)[names(Workflow$.__enclos_env__$private$biasFieldsSpecify) %in% names(speciesDataset)]

      for (biasName in specifySubset) {

        initializeModel$spatialFields$biasFields[[biasName]] <- Workflow$.__enclos_env__$private$biasFieldsSpecify[[biasName]]

      }

    }

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
                                                                mask = .__mask.__,
                                                                dims = predictionDim),
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

  if ('Bias' %in% Oputs) {

    if (biasIn) {

    if (!Quiet) message('\nProducing bias predictions:\n\n')
    .__mask.__ <- as(Workflow$.__enclos_env__$private$Area, 'Spatial')
    biasPreds <- predict(PSDMsMOdel,
                         data = inlabru::fm_pixels(mesh = .__mesh.__,
                         mask = .__mask.__,
                         dims = predictionDim),
                         biasfield = TRUE)

    if (saveObjects) {

      if (!Quiet)  message('\nSaving predictions object:', '\n\n')
      saveRDS(object = biasPreds, file = paste0(modDirectory,'/', speciesNameInd, '/biasPreds.rds'))

    } else outputList[[speciesNameInd]][['Bias']] <- biasPreds

    }

  }

  }

  if ('Richness' %in% Oputs) {

    message('Richness output not available yet!')

    spData <- append(Workflow$.__enclos_env__$private$dataGBIF,
                     Workflow$.__enclos_env__$private$dataStructured) #Check

    richSetup <- PointedSDMs::intModel(speciesDataset, Mesh = .__mesh.__, Projection = .__proj.__, Coordinates = .__coordinates.__,
                          responsePA = .__responsePA.__, responseCounts = .__responseCounts.__,
                          trialsPA = .__trialsName.__,
                          pointsIntercept = .__pointsIntercept.__ ,
                          copyModel = .__copyModel.__, speciesName = workflow$.__enclos_env__$private$speciesName,
                          speciesSpatial = 'shared', ##WHICH ONE??
                          pointsSpatial = 'copy', speciesIndependent = TRUE,
                          speciesEffects = list(Intercept = FALSE, Environmental = TRUE),
                          spatialCovariates = spatCovs)

    if (!is.null(Workflow$.__enclos_env__$private$biasNames)) {

      richSetup$addBias(datasetNames = Workflow$.__enclos_env__$private$biasNames, copyModel = Workflow$.__enclos_env__$private$biasFieldsCopy)

      if (!is.null(Workflow$.__enclos_env__$private$biasFieldsSpecify)) {

        for (biasName in names(Workflow$.__enclos_env__$private$biasFieldsSpecify)) {

          richSetup$spatialFields$biasFields[[biasName]] <- Workflow$.__enclos_env__$private$biasFieldsSpecify[[biasName]]

        }


      }


      richModel <- PointedSDMs::fitISDM(data = richSetup,
                                        options = Workflow$.__enclos_env__$private$optionsINLA)

      ##Predict and plot


    }

  }

  if (!saveObjects) return(outputList)


  }
