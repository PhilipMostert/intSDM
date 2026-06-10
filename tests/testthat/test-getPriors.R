testthat::test_that('getPriors obtains the correct priors from an INLA object', {

  skip_on_cran()

  inlaOptions <- list(control.inla = list(int.strategy = 'eb', diagonal = 1))

  test_data <- readRDS(system.file('extdata/test_data.rds', package = 'intSDM'))

  POpoints <- test_data$POpoints
  PApoints <- test_data$PApoints
  Mesh <- test_data$Mesh

  proj <- '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=km +no_defs'

  countries <- st_as_sf(geodata::world(path = tempdir()))

  countries <- countries[countries$NAME_0 %in% c('Norway'),]
  countries <- st_transform(countries, proj)
  species <- c('Fraxinus excelsior')
  workflow <- try(startWorkflow(Species = species,
                                saveOptions = list(projectName = 'testthatexample', projectDirectory = './'),
                                Projection = proj, Countries = 'Norway',
                                Quiet = TRUE, Save = FALSE))

  workflow$addStructured(dataStructured = POpoints, datasetType = 'PO', datasetName = 'PO', speciesName = 'name')
  workflow$addStructured(dataStructured = PApoints, datasetType = 'PA', datasetName = 'PA', speciesName = 'name', responseName = 'pres')

  workflow$workflowOutput('Model')
  workflow$addMesh(Mesh) #200000

  workflow$modelOptions(ISDM = list(pointsSpatial = 'shared'))

  workflow$specifyPriors(effectNames = 'Intercept',
                         Mean = 0,
                         Precision = 1)

  workflow$specifySpatial(prior.range = c(100, 0.1), prior.sigma = c(1, 0.2))

  model <- sdmWorkflow(Workflow = workflow, inlaOptions = inlaOptions)

  priors <- getPriors(model$Fraxinus_excelsior$Model)

  expect_setequal(names(priors), c('fixedEffects', 'randomEffects'))

  expect_setequal(row.names(priors$fixedEffects), c('PO_intercept', 'PA_intercept'))

  expect_all_equal(priors$fixedEffects[,1], 0)
  expect_all_equal(priors$fixedEffects[,2], 1)

  expect_setequal(row.names(priors$fixedEffects), c('PO_intercept', 'PA_intercept'))

  expect_equal(priors$randomEffects$shared_spatial$prior, 'pcmatern')
  expect_equal(priors$randomEffects$shared_spatial$values$Range, c(range = 100, prob = 0.1))
  expect_equal(priors$randomEffects$shared_spatial$values$StDev, c(sigma = 1, prob = 0.2))



}
  )
