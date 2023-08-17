##First set up workflow
library(lwgeom)
proj <- '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
species <- 'Fraxinus excelsior'
workflow <- startWorkflow(Species = species,
                          saveOptions = list(projectName = 'testthatexample'),
                          Projection = proj,
                          Quiet = TRUE, Save = FALSE)

testthat::test_that('Test that addArea correctly adds the correct area to the model', {
  skip_on_cran()

  #Check that countryName works
  expect_error(workflow$addArea(), 'One of object or countryName is required.')

  workflow$addArea(countryName = c('Sweden', 'Norway'))
  expect_setequal(class(workflow$.__enclos_env__$private$Area), c('sf', 'data.frame'))
  expect_setequal(workflow$.__enclos_env__$private$Area$NAME_ENGL, c('Sweden', 'Norway'))
  expect_identical(st_crs(workflow$.__enclos_env__$private$Area)[2], st_crs(proj)[2])

  ##Check adding an area as an object
  workflow2 <- startWorkflow(Species = species,
                            saveOptions = list(projectName = 'testthatexample'),
                            Projection = proj,
                            Quiet = TRUE, Save = FALSE)
  #Obtain object
  countries <<- giscoR::gisco_countries[giscoR::gisco_countries$NAME_ENGL %in% c('Sweden', 'Norway'), ]
  workflow2$addArea(Object = countries)
  #expect_identical(workflow2$.__enclos_env__$private$Area, countries) Won't be identical due to change in CRS

  countriesSP <- as(countries, 'Spatial')
  workflow2$addArea(Object = countriesSP)

  expect_error(workflow2$addArea(Object = data.frame(x = runif(100), y = runif(100))), 'Object needs to be a sp or sf object.')


})

testthat::test_that('Test that addGBIF correctly adds the correct data to the model', {

  skip_on_cran()

  expect_error(workflow$addGBIF(), 'Please provide a name to give your dataset using datasetName.')
  expect_error(workflow$addGBIF(Species = 'Not_provided', datasetName = 'TEST'), 'Species provided not specified in startWorkflow().')

  workflow$addGBIF(datasetName = 'GBIFTEST')

  expect_equal(workflow$.__enclos_env__$private$classGBIF$Fraxinus_excelsior$GBIFTEST, 'PO')
  expect_equal(class(workflow$.__enclos_env__$private$dataGBIF), 'list')
  expect_equal(names(workflow$.__enclos_env__$private$dataGBIF), 'Fraxinus_excelsior')
  expect_equal(names(workflow$.__enclos_env__$private$dataGBIF$Fraxinus_excelsior), 'GBIFTEST')
  expect_setequal(class(workflow$.__enclos_env__$private$dataGBIF$Fraxinus_excelsior$GBIFTEST), c('sf', 'data.frame'))

  #Change dataset type to PA
  workflow$addGBIF(datasetName = 'GBIFTEST2', datasetType = 'PA')

  expect_equal(workflow$.__enclos_env__$private$classGBIF$Fraxinus_excelsior$GBIFTEST2, 'PA')
  expect_setequal(class(workflow$.__enclos_env__$private$dataGBIF$Fraxinus_excelsior$GBIFTEST2), c('sf', 'data.frame'))
  expect_true(unique(workflow$.__enclos_env__$private$dataGBIF$Fraxinus_excelsior$GBIFTEST2$occurrenceStatus) %in% c(0,1))

  ##Change to Counts and check that NAs are removed.
  expect_warning(workflow$addGBIF(datasetName = 'GBIFTEST3', datasetType = 'Counts', limit = 1000), 'Removing reccords with NA individualCount values')
  expect_setequal(class(workflow$.__enclos_env__$private$dataGBIF$Fraxinus_excelsior$GBIFTEST3), c('sf', 'data.frame'))
  expect_true(sum(is.na(workflow$.__enclos_env__$private$dataGBIF$GBIFTEST3$individualCount)) == 0)

  ##Check assign2global
  workflow$addGBIF(datasetName = 'GBIFTEST4', assign2Global = TRUE)

  expect_true('GBIFTEST4' %in% ls(envir = globalenv()))
  #rm(GBIFTEST4, envir = globalenv())


})

testthat::test_that('Test that addCovariate correctly adds the desired covariate to the model', {

  skip_on_cran()

  skip(message = 'geodata is down for now')
  expect_error(workflow$addCovariates(), 'One of object or worldClim is required.')

  covariateWorkflow <- startWorkflow(Species = species,
                                     saveOptions = list(projectName = 'testthatexample',
                                                        projectDirectory = './tests'),
                                     Projection = proj,
                                     Quiet = TRUE, Save = TRUE)

  covariateWorkflow$addArea(countryName = c('Norway', 'Sweden'))

  expect_error(covariateWorkflow$addCovariates(worldClim = 'prec', Months = c('Monday', 'Tuesday')), 'Month provided is not valid.')
  expect_error(covariateWorkflow$addCovariates(worldClim = 'depth', Months = c('June', 'July', 'August')), 'worldClim argument is not a valid option.')
  expect_error(covariateWorkflow$addCovariates(worldClim = c('prec', 'bio'), Months = c('June', 'July', 'August')), 'Please only add one worldClim variable at a time.')

  covariateWorkflow$addCovariates(worldClim = 'prec', Months = c('June', 'July', 'August'))

  expect_equal(names(covariateWorkflow$.__enclos_env__$private$Covariates), 'prec')
  expect_equal(class(covariateWorkflow$.__enclos_env__$private$Covariates$prec)[1], 'SpatRaster')
  expect_true(length(covariateWorkflow$.__enclos_env__$private$Covariates) == 1)
  expect_true(names(covariateWorkflow$.__enclos_env__$private$Covariates$prec) == 'prec')

  unlink('./tests/testthatexample', recursive = TRUE)

  ##Test adding own covariate layer
  NorSwe <- geodata::worldclim_country(country = c('Norway', 'Sweden'), var = 'tavg', path = './tests/testthatremove')
  Port <- geodata::worldclim_country(country = c('Portugal'), var = 'tavg', path = './tests/testthatremove')

  expect_error(covariateWorkflow$addCovariates(Object = Port[[1]]), 'The covariate provided and the area specified do not match.')
  expect_error(covariateWorkflow$addCovariates(Object = NorSwe), 'Please provide each covariate into the workflow as their own object.')


  covariateWorkflow$addCovariates(Object = NorSwe[[1]])
  expect_setequal(names(covariateWorkflow$.__enclos_env__$private$Covariates), c("prec", "NOR_wc2.1_30s_tavg_1"))

  NorSwe <- as(NorSwe, 'Raster')
  covariateWorkflow$addCovariates(Object = NorSwe[[1]])
  expect_equal((class(covariateWorkflow$.__enclos_env__$private$Covariates$NOR_wc2.1_30s_tavg_1))[1], 'SpatRaster')

  NorSwe <- as(NorSwe, 'SpatialPixelsDataFrame')
  covariateWorkflow$addCovariates(Object = NorSwe[, 1])
  expect_equal((class(covariateWorkflow$.__enclos_env__$private$Covariates$NOR_wc2.1_30s_tavg_1))[1], 'SpatRaster')

  unlink('./tests/testthatremove', recursive = TRUE)


})

testthat::test_that('addStructured can add the data correctly to the model', {

  skip_on_cran()


  workflow <- startWorkflow(Species = species,
                            saveOptions = list(projectName = 'testthatexample'),
                            Projection = proj,
                            Quiet = TRUE, Save = FALSE)
  workflow$addArea(countryName = c('Sweden', 'Norway'))

  #Simulate PA datasets.
  dataPA <- st_as_sf(st_sample(x = countries, size = 100))
  st_geometry(dataPA) <- 'geometry'
  dataPA$Presence <- sample(c(0,1), 100, replace = TRUE)
  dataPA$species <- sample(c('Bird', 'Fish'), 100, replace = TRUE)

  expect_error(workflow$addStructured(), 'dataStructured needs to be provided')
  expect_error(workflow$addStructured(dataStructured = dataPA), 'speciesName cannot be missing.')
  expect_error(workflow$addStructured(dataStructured = dataPA, speciesName = 'species'), 'datasetType needs to be one of "PO", "PA" or "Counts".')
  expect_error(workflow$addStructured(dataStructured = dataPA, datasetType = 'PA', speciesName = 'species'), 'responseName cannot be missing for PA and counts datasets.')

  expect_error(workflow$addStructured(dataStructured = dataPA, datasetType = 'PA', responseName = 'Presence', speciesName = 'species'))

  dataPA$species <- sample(c('Bird', species), 100, replace = TRUE)
  expect_warning(workflow$addStructured(dataStructured = dataPA, datasetType = 'PA', responseName = 'Presence', speciesName = 'species'))

  expect_equal(names(workflow$.__enclos_env__$private$dataStructured), sub(" ", "_", species))
  expect_equal(names(workflow$.__enclos_env__$private$dataStructured$Fraxinus_excelsior), 'dataPA')
  expect_setequal(names(workflow$.__enclos_env__$private$dataStructured$Fraxinus_excelsior$dataPA), c('geometry', 'occurrenceStatus', 'speciesName'))

  #Simulate repeated PA

  dataPA2 <- st_as_sf(st_sample(x = countries, size = 100))
  st_geometry(dataPA2) <- 'geometry'
  dataPA2$Presence <- sample(c(0,1), 100, replace = TRUE)
  dataPA2$trials <- sample(c(2,3), 100, replace = TRUE)
  dataPA2$species <- species

  workflow$addStructured(dataStructured = dataPA2, datasetType = 'PA', speciesName = 'species',
                         responseName = 'Presence', trialsName = 'trials')

  expect_setequal(names(workflow$.__enclos_env__$private$dataStructured$Fraxinus_excelsior), c('dataPA', 'dataPA2'))
  expect_setequal(names(workflow$.__enclos_env__$private$dataStructured$Fraxinus_excelsior$dataPA2), c('geometry', 'occurrenceStatus', 'numTrials', 'speciesName'))

  #Simulate Counts data
  dataCounts <- st_as_sf(st_sample(x = countries, size = 100))
  st_geometry(dataCounts) <- 'geometry'
  dataCounts$num <- rpois(n = 100, lambda = 10)
  dataCounts$species <- species

  workflow$addStructured(dataStructured = dataCounts, datasetType = 'Counts', responseName = 'num', speciesName = 'species')

  expect_setequal(names(workflow$.__enclos_env__$private$dataStructured$Fraxinus_excelsior), c('dataPA', 'dataPA2', 'dataCounts'))
  expect_setequal(names(workflow$.__enclos_env__$private$dataStructured$Fraxinus_excelsior$dataCounts), c('geometry', 'individualCount', 'speciesName'))

  #Add an sp dataset
  #dataPASP <- as(dataPA, 'Spatial')
  #workflow$addStructured(dataStructured = dataPASP, datasetType = 'PA', responseName = 'Presence', speciesName = 'species')
  #expect_setequal(class(workflow$.__enclos_env__$private$dataStructured$Fraxinus_excelsior$dataPASP), c('sf', 'data.frame'))

  #dataSP <- as(st_as_sf(st_sample(x = countries, size = 100)), 'Spatial')
  #dataSP$species <- species
  #workflow$addStructured(dataStructured = dataSP, datasetType = 'PO', speciesName = 'species')
  #expect_setequal(names(workflow$.__enclos_env__$private$dataStructured$Fraxinus_excelsior$dataSP), c('geometry', 'speciesName'))
  #expect_setequal(class(workflow$.__enclos_env__$private$dataStructured$Fraxinus_excelsior$dataSP), c('sf', 'data.frame'))

  #Add a data.frame object
  dataFrame <- st_transform(dataPA, proj)
  dataFrame$Presence <- NULL
  dataFrame <- data.frame(st_coordinates(dataFrame))
  dataFrame$species <- species

  workflow$addStructured(dataStructured = dataFrame, datasetType = 'PO', coordinateNames = c('X', 'Y'), speciesName = 'species')
  expect_setequal(class(workflow$.__enclos_env__$private$dataStructured$Fraxinus_excelsior$dataFrame), c('sf', 'data.frame'))

  #Add data not in boundary
  dataNotIn <- st_as_sf(st_sample(x = giscoR::gisco_countries[giscoR::gisco_countries$NAME_ENGL == 'Portugal',], size = 100))
  dataNotIn$species <- species
  expect_error(workflow$addStructured(dataStructured = dataNotIn, datasetType = 'PO', speciesName = 'species'), 'Dataset provided has no reccords over the boundary.')

  workflow <<- workflow

  unlink('./tests/testthatexample', recursive = TRUE)

})

testthat::test_that('addMesh correctly adds the mesh to the model', {

  skip_on_cran()

  expect_error(workflow$addMesh())

  meshObject <<- INLA::inla.mesh.2d(boundary = inlabru::fm_sp2segment(countries),
                             max.edge = 200000,
                             offset = 10,
                             crs = proj)

  workflow$addMesh(Object = meshObject)
  expect_equal(class(workflow$.__enclos_env__$private$Mesh), 'inla.mesh')

  workflow$addMesh(max.edge = 200000, cutoff = 3)
  expect_equal(class(workflow$.__enclos_env__$private$Mesh), 'inla.mesh')

})

testthat::test_that('crossValidation correctly specifies the correct cross-validation method', {

  expect_error(workflow$crossValidation())

  expect_error(workflow$crossValidation(Method = 'kfold'), 'Output needs to be at least one of: spatialBlock, Loo.')

  workflow$crossValidation(Method = c('spatialBlock', 'Loo'))

  expect_setequal(workflow$.__enclos_env__$private$CVMethod, c("spatialBlock", "Loo"))

  expect_error(workflow$crossValidation(Method = 'spatialBlock', blockOptions = list(k =2)), 'Please provide both k and rows_cols in blockOptions.')

  workflow$crossValidation(Method = 'spatialBlock', blockOptions = list(k =2, rows_cols = c(4,5)))

  expect_equal(class(workflow$.__enclos_env__$private$blockOptions), 'list')

  expect_equal(workflow$.__enclos_env__$private$blockOptions$k, 2)

  expect_equal(workflow$.__enclos_env__$private$blockOptions$rows_cols, c(4,5))


})

testthat::test_that('modelOptions correctly adds options', {

  skip_on_cran()

  expect_error(workflow$modelOptions(INLA = 'control.strategy = list()'), 'INLA needs to be a list of INLA arguments to specify the model.')

  expect_error(workflow$modelOptions(ISDM = list(marks = TRUE)), 'ISDM needs to be a named list with at least one of the following options: "pointCovariates", "pointsIntercept", "pointsSpatial" or "copyModel".')
  expect_error(workflow$modelOptions(ISDM = list(pointsSpatial = FALSE, marks = TRUE)), 'ISDM needs to be a named list with at least one of the following options: "pointCovariates", "pointsIntercept", "pointsSpatial" or "copyModel".')

  workflow$modelOptions(ISDM = list(pointsSpatial = 'copy', copyModel = list(beta = list(fixed = TRUE))))
  workflow$modelOptions(INLA = list(control.inla = list(int.strategy = 'eb')))

  expect_setequal(names(workflow$.__enclos_env__$private$optionsISDM), c('pointsSpatial', 'copyModel'))

  expect_true(workflow$.__enclos_env__$private$optionsISDM$copyModel$beta$fixed == TRUE)
  expect_true(workflow$.__enclos_env__$private$optionsISDM$pointsSpatial == 'copy')

  expect_true(workflow$.__enclos_env__$private$optionsINLA$control.inla$int.strategy == 'eb')

})

testthat::test_that('specifySpatial correctly specifies the spatial fields', {

  skip_on_cran()

  expect_error(workflow$specifySpatial(), 'Please provide arguments to customize the INLA spde object using the ... argument.')
  workflowNoMesh <- startWorkflow(Species = species,
                                  saveOptions = list(projectName = 'testthatexample'),
                                  Projection = proj,
                                  Quiet = TRUE, Save = FALSE)

  expect_error(workflowNoMesh$specifySpatial(prior.range = c(1,0.1), prior.sigma = c(1, 0.2)))

  workflow$specifySpatial(prior.range = c(1, 0.1), prior.sigma = c(1, 0.2))

  expect_setequal(class(workflow$.__enclos_env__$private$sharedField), c("inla.spde2", "inla.spde", "inla.model.class"))

  })

testthat::test_that('biasFields correctly adds the bias field', {

  skip_on_cran()

  expect_error(workflow$biasFields(), 'argument "datasetName" is missing, with no default')

  expect_error(workflow$biasFields(datasetName = 'NotIn'), 'Dataset specified for bias field not included in the workflow.')

  workflow$biasFields(datasetName = c('dataCounts', "dataFrame"))
  expect_setequal(workflow$.__enclos_env__$private$biasNames, c("dataCounts", "dataFrame"))

  workflow$biasFields(datasetName = c('dataPA2', 'dataPA'))

  expect_setequal(workflow$.__enclos_env__$private$biasNames, c("dataCounts", "dataFrame", 'dataPA2', 'dataPA'))

  workflow$biasFields(datasetName = c('dataPA2', 'dataPA'), prior.range = c(1, 0.1), prior.sigma = c(1, 0.2))

  expect_setequal(names(workflow$.__enclos_env__$private$biasFieldsSpecify), c('dataPA2', 'dataPA'))

  expect_setequal(class(workflow$.__enclos_env__$private$biasFieldsSpecify$dataPA2), c("inla.spde2", "inla.spde", "inla.model.class"))
  expect_setequal(class(workflow$.__enclos_env__$private$biasFieldsSpecify$dataPA), c("inla.spde2", "inla.spde", "inla.model.class"))



})

testthat::test_that('workflowOutput gives the correct output', {

  skip_on_cran()

  expect_error(workflow$workflowOutput(), 'argument "Output" is missing, with no default')

  expect_error(workflow$workflowOutput('Trendline'), 'Output needs to be at least one of: Model, Predictions, Maps or Cross-validation.')

  workflow$workflowOutput(c('Model', 'Maps', 'Cross-validation'))

  })



