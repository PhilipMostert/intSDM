sdmWorkflow <- function(workflow = NULL,
                        saveObjects = NULL,
                        Report = FALSE) {

  if (is.null(workFlow$.__enclos_env__$private$Mesh)) stop('An inla.mesh object is required before any analysis is completed. Please add using the `.$addMesh` function.')
  if (is.null(workflow$.__enclos_env__$private$Output)) stop('A model output needs to be specified before any analysis is completed. Please add using the `.$workflowOutput` function.')

  Oputs <- workflow$.__enclos_env__$private$Output

  if (is.null(workflow$.__enclos_env__$private$CVMethod) && 'Cross-Validation' %in% workflow$.__enclos_env__$private$Output) stop('Cross-validation specified as model output but no method provided. Please specify cross-validation method using the `.$crossValidation` function.')

  dataUsed <- append(workflow$.__enclos_env__$private$dataGBIF,
                     workflow$.__enclos_env__$private$dataStructured)

  if (length(dataUsed) == 0) stop('No data added to the model. Please add data using `.$addGBIF` or `.$addStructured`.')

  if (length(workFlow$optionsISDM) > 1) {

    ##Add options here


  } else {

    ##Add default options here

  }

  initializeModel <- PointedSDMs::intModel()

  if ('Cross-validation' %in% Oputs && 'spatialBlock' %in%workflow$.__enclos_env__$private$CVMethod) {

    ##Somehow add blockCV model nicely

  }

  message('Estimating ISDM:\n\n')

  PSDMsMOdel <- PointedSDMs::fitISDM(initializeModel,
                                     options = workflow$.__enclos_env$private$optionsINLA)

  if ('Model' %in% Oputs) {

    message('Saving Model object to:', '\n\n')
    saveRDS(object = PSDMsMOdel)


  }

  if ('Cross-validation' %in% Oputs) {


    if ('spatialBlock' %in %workflow$.__enclos_env__$private$CVMethod) {

      message('Completing spatial block cross-validation:\n\n')
      spatialBlockCV <- PointedSDMs::blockedCV(initializeModel)
      message('Saving spatial blocked cross-validation object to:', '\n\n')
      saveRDS(object = spatialBlockCV); rm(spatialBlockCV)

    }

    if ('Loo' %in% workflow$.__enclos_env__$private$CVMethod) {

      message('Completing leave-one-out cross-validation:\n\n')
      LooCV <- PointedSDMs::datasetOut()
      message('Saving leave-one-out cross-validation object to:', '\n\n')
      saveRDS(object = LooCV); rm(LooCV)

    }

    if ('Predictions' %in% Oputs) {

      message('Predicting model:\n\n')
      Predictions <- predict(model, data = inlabru::pixels(mesh, mask), predictor = TRUE)
      message('Saving predictions object to:', '\n\n')
      saveRDS(object = Predictions)

    }

    if ('Maps' %in% Oputs) {

      message("Plotting Maps:\n\n")
      plot(Predictions)
      #Save object here

    }

    ##Some summary output here
     #Should we delete all the objects after use?

  }

  }
