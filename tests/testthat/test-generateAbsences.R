testthat::test_that('generateAbsences correctly creates absences for the data.', {

  ##First set up workflow
  skip_on_cran()

  proj <- '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
  species <- c('Fraxinus excelsior', 'Ulmus glabra', 'Arnica montana')
  workflow <- startWorkflow(Species = species,
                            saveOptions = list(projectName = 'testthatexample'),
                            Projection = proj,
                            Quiet = TRUE, Save = FALSE)
  workflow$addArea(countryName = c('Sweden', 'Norway'))

  workflow$addGBIF(datasetType = 'PO', limit = 50, datasetName = 'PO')

  workflow$addGBIF(datasetType = 'PA', datasetName = 'PA', generateAbsences = FALSE)

  paData <- lapply(workflow$.__enclos_env__$private$dataGBIF, function(x) x[['PA']])

  workflow$addGBIF(datasetType = 'PA', datasetName = 'PA', generateAbsences = TRUE)

  expect_true(all(names(workflow$.__enclos_env__$private$dataGBIF) %in% sub(" ", '_', species)))

  expect_true(all(unlist(lapply(workflow$.__enclos_env__$private$dataGBIF, function(x) names(x))) %in% c('PO', 'PA')))

  expect_true(nrow(paData$Fraxinus_excelsior) < nrow(workflow$.__enclos_env__$private$dataGBIF$Fraxinus_excelsior$PA))
  expect_true(nrow(paData$Ulmus_glabra) < nrow(workflow$.__enclos_env__$private$dataGBIF$Ulmus_glabra$PA))
  expect_true(nrow(paData$Arnica_montana) < nrow(workflow$.__enclos_env__$private$dataGBIF$Arnica_montana$PA))

  expect_true(all(is.na(data.frame(workflow$.__enclos_env__$private$dataGBIF$Fraxinus_excelsior$PA)[(nrow(paData$Fraxinus_excelsior) +1: nrow(workflow$.__enclos_env__$private$dataGBIF$Arnica_montana$PA)),
  !names(workflow$.__enclos_env__$private$dataGBIF$Fraxinus_excelsior$PA) %in% c('species', 'networkKeys','occurrenceStatus', 'geometry')])))

  expect_true(all(is.na(data.frame(workflow$.__enclos_env__$private$dataGBIF$Ulmus_glabra$PA)[(nrow(paData$Ulmus_glabra) +1: nrow(workflow$.__enclos_env__$private$dataGBIF$Ulmus_glabra$PA)),
                                                                                                    !names(workflow$.__enclos_env__$private$dataGBIF$Ulmus_glabra$PA) %in% c('species', 'networkKeys','occurrenceStatus', 'geometry')])))

  expect_true(all(is.na(data.frame(workflow$.__enclos_env__$private$dataGBIF$Arnica_montana$PA)[(nrow(paData$Arnica_montana) +1: nrow(workflow$.__enclos_env__$private$dataGBIF$Arnica_montana$PA)),
                                                                                              !names(workflow$.__enclos_env__$private$dataGBIF$Arnica_montana$PA) %in% c('species', 'networkKeys','occurrenceStatus', 'geometry')])))




})
