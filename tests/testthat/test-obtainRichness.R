testthat::test_that('obtainRichness can produce an sf object of species richness', {

  skip_on_cran()
  skip_if_not_installed("R.utils")
  skip_if_not(local_testthat_geodata_path())

  library(R.utils)

  inlaOptions <- list(control.inla = list(int.strategy = 'eb', diagonal = 1))

  test_data <- readRDS(system.file('extdata/test_data.rds', package = 'intSDM'))

  POpoints <- test_data$POpoints
  PApoints <- test_data$PApoints
  Mesh <- test_data$Mesh

  proj <- '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=km +no_defs'
  countries <- st_as_sf(geodata::world(path = geodata::geodata_path()))
  countries <- countries[countries$NAME_0 %in% c('Norway'),]
  countries <- st_transform(countries, proj)
  species <- c('Fraxinus excelsior')
  projDir <- tempfile("intSDM_obtainRichness")
  on.exit(unlink(projDir, recursive = TRUE))

  workflow <- try(startWorkflow(Species = species,
                                saveOptions = list(projectName = 'testthatexample', projectDirectory = projDir),
                                Projection = proj, Countries = 'Norway', Richness = TRUE,
                                Quiet = TRUE, Save = FALSE))

  if (inherits(workflow, 'try-error')) {


    workflow <- startWorkflow(Species = species,
                              saveOptions = list(projectName = 'testthatexample', projectDirectory = projDir),
                              Projection = proj, Quiet = TRUE, Save = FALSE, Richness = TRUE)


    workflow$addArea(Object = countries)

  }

  workflow$addStructured(dataStructured = POpoints, datasetType = 'PO', datasetName = 'PO', speciesName = 'name')
  workflow$addStructured(dataStructured = PApoints, datasetType = 'PA', datasetName = 'PA', speciesName = 'name', responseName = 'pres')

  workflow$workflowOutput(c('Model'))

  covs <- terra::rast(system.file('extdata/vignette_covariates.tif', package = 'intSDM'))

  workflow$addCovariates(covs)

  workflow$addMesh(Object = Mesh) #200000
  workflow$modelOptions(Richness = list(predictionIntercept = 'PA'))
  #Make data if NA
  model <- sdmWorkflow(Workflow = workflow, inlaOptions = inlaOptions)

  ##Try wrong modelObject
  Rich <- expect_error(obtainRichness(modelObject = model), 'modelObject needs to be a modSpecies object obtained from the PointedSDMs function fitISDM.')
  Rich <- expect_error(obtainRichness(modelObject = model$RichnessModel,
                                      predictionData = fmesher::fm_pixels(workflow$.__enclos_env__$private$Mesh)),'predictionIntercept cannot be missing.')
  Rich <- expect_error(obtainRichness(modelObject = model$RichnessModel,
                                      predictionData = fmesher::fm_pixels(workflow$.__enclos_env__$private$Mesh),
                                      predictionIntercept = 'wrong'),'predictionIntercept needs to be the name of a dataset included in modelObject.')

  predDat <- fmesher::fm_pixels(workflow$.__enclos_env__$private$Mesh)
  predDat$Cov1 <- rnorm(nrow(predDat))
  predDat$Cov2 <- rnorm(nrow(predDat))
  Rich <- obtainRichness(modelObject = model$RichnessModel,
                       predictionData = predDat,
                       predictionIntercept = 'PA')

  expect_identical(class(Rich), 'list')
  expect_setequal(names(Rich), c('Richness', 'Probabilities'))
  expect_setequal(names(Rich$Probabilities), gsub(' ', '_', species))

  })
