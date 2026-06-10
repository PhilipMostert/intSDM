testthat::test_that('sdmWorkflow produces the correct output given different Workflow situations.', {

  skip_on_cran()
  skip_if_not(local_testthat_geodata_path())

  projectDir <- tempfile("intSDM_sdmWorkflow_test_")
  on.exit(unlink(projectDir, recursive = TRUE))
  dir.create(projectDir, recursive = TRUE, showWarnings = FALSE)

  ##Create different workflows here:
  #1. Just GBIF data
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
  workflow <- try(startWorkflow(Species = species,
                                saveOptions = list(projectName = 'testthatexample', projectDirectory = projectDir),
                                Projection = proj, Countries = 'Norway',
                                Quiet = TRUE, Save = TRUE))

  if (inherits(workflow, 'try-error')) {


    workflow <- startWorkflow(Species = species,
                              saveOptions = list(projectName = 'testthatexample', projectDirectory = projectDir),
                              Projection = proj, Quiet = TRUE, Save = TRUE)


    workflow$addArea(Object = countries)

  }

  workflow$addStructured(dataStructured = POpoints, datasetType = 'PO', datasetName = 'PO', speciesName = 'name')
  workflow$addStructured(dataStructured = PApoints, datasetType = 'PA', datasetName = 'PA', speciesName = 'name', responseName = 'pres')

  expect_error(sdmWorkflow(Workflow = workflow)) #Test no output given
  workflow$workflowOutput('Model')
  expect_error(sdmWorkflow(Workflow = workflow)) #Test no mesh provided
  workflow$addMesh(Mesh)

  workflow$modelOptions(ISDM = list(pointsSpatial = 'shared'))
  # May generate a warning (under rstudio);
  # 'package:stats' may not be available when loading
  suppressWarnings(sdmWorkflow(Workflow = workflow, inlaOptions = inlaOptions))
  expect_true(all(c(dir.exists(file.path(projectDir, 'testthatexample', 'Fraxinus_excelsior')))))

  expect_true(all(c(file.exists(file.path(projectDir, 'testthatexample', 'Fraxinus_excelsior' , 'intModel.rds')))))

  Fraxinus_excelsior <- readRDS(file = file.path(projectDir, 'testthatexample', 'Fraxinus_excelsior', 'intModel.rds'))
  expect_setequal(rownames(Fraxinus_excelsior$summary.fixed), c("PO_intercept", "PA_intercept"))
  expect_equal(as.character(Fraxinus_excelsior$componentsJoint)[2],
               "-1 + shared_spatial(main = geometry, model = shared_field) + PO_intercept(1) + PA_intercept(1)")
  rm(Fraxinus_excelsior)
  biasCopyonCRAN <- FALSE
  if (biasCopyonCRAN) {

    biasWorkflow <- startWorkflow(Species = species,
                                  saveOptions = list(projectName = 'testthatexample2', projectDirectory = projectDir),
                                  Projection = proj,
                                  Quiet = TRUE, Save = FALSE)

    biasWorkflow$addArea(Object = countries)


    biasWorkflow$addStructured(dataStructured = POpoints, datasetType = 'PO', datasetName = 'PO', speciesName = 'name')
    biasWorkflow$addStructured(dataStructured = PApoints, datasetType = 'PA', datasetName = 'PA', speciesName = 'name', responseName = 'pres')

    biasWorkflow$workflowOutput('Model')
    biasWorkflow$addMesh(Mesh) #200000
    biasWorkflow$biasFields('PO')
    biasWorkflow$modelOptions(ISDM = list(pointsSpatial = 'shared'))
    biasMod <- sdmWorkflow(biasWorkflow, inlaOptions = inlaOptions)

    expect_setequal(names(biasMod$Fraxinus_excelsior$Model$summary.random), c("shared_spatial", "PO_biasField"))
    expect_setequal(class(biasMod$Fraxinus_excelsior$Model), c("modISDM", "bru", "iinla", "inla"))
    rm(biasWorkflow)

  }
  copyWorkflow <- startWorkflow(Species = species,
                                saveOptions = list(projectName = 'testthatexample3', projectDirectory = projectDir),
                                Projection = proj,
                                Quiet = TRUE, Save = FALSE)

  copyWorkflow$addArea(Object = countries)

  copyWorkflow$addStructured(dataStructured = POpoints, datasetType = 'PO', datasetName = 'PO', speciesName = 'name')
  copyWorkflow$addStructured(dataStructured = PApoints, datasetType = 'PA', datasetName = 'PA', speciesName = 'name', responseName = 'pres')

  copyWorkflow$workflowOutput('Model')
  copyWorkflow$addMesh(Mesh) #200000
  copyWorkflow$modelOptions(ISDM = list(pointsSpatial = 'copy'))
  copyWorkflow$specifyPriors(copyModel = list(beta = list(fixed = TRUE)))
  copyMod <- sdmWorkflow(copyWorkflow, inlaOptions = inlaOptions)

  expect_setequal(names(copyMod$Fraxinus_excelsior$Model$summary.random), c("PO_spatial", "PA_spatial"))
  expect_equal(as.character(copyMod$Fraxinus_excelsior$Model$componentsJoint)[2],
               "-1 + PO_spatial(main = geometry, model = PO_field) + PA_spatial(main = geometry, copy = \"PO_spatial\", hyper = list(beta = list(fixed = TRUE))) + PO_intercept(1) + PA_intercept(1)")

  ##Test Richness model
  species <- c('Fraxinus excelsior')
  workflow <- try(startWorkflow(Species = species,
                                saveOptions = list(projectName = 'richness', projectDirectory = projectDir),
                                Projection = proj, Countries = 'Norway', Richness = TRUE,
                                Quiet = TRUE, Save = TRUE))

  if (inherits(workflow, 'try-error')) {


    workflow <- startWorkflow(Species = species,
                              saveOptions = list(projectName = 'richness', projectDirectory = projectDir),
                              Projection = proj, Quiet = TRUE, Save = TRUE, Richness = TRUE)


    workflow$addArea(Object = countries)

  }

  workflow$addStructured(dataStructured = POpoints, datasetType = 'PO', datasetName = 'PO', speciesName = 'name')
  workflow$addStructured(dataStructured = PApoints, datasetType = 'PA', datasetName = 'PA', speciesName = 'name', responseName = 'pres')

  expect_error(sdmWorkflow(Workflow = workflow)) #Test no output given

  workflow$workflowOutput(c('Model', 'Predictions'))

  expect_error(sdmWorkflow(Workflow = workflow)) #Test no mesh provided

  workflow$addMesh(Mesh) #200000

  workflow$modelOptions(ISDM = list(pointsSpatial = 'shared'))
  expect_error(sdmWorkflow(Workflow = workflow))
  workflow$modelOptions(Richness = list(predictionIntercept = 'PA'))
  # May generate a warning (under rstudio);
  # 'package:stats' may not be available when loading
  suppressWarnings(sdmWorkflow(Workflow = workflow, inlaOptions = inlaOptions))

  expect_true(all(c(file.exists(file.path(projectDir, 'richness', 'richnessModel.rds')))))

  expect_true(all(c(file.exists(file.path(projectDir, 'richness', 'richnessPredictions.rds')))))

  RichModel <- readRDS(file = file.path(projectDir, 'richness', 'richnessModel.rds'))
  expect_setequal(rownames(RichModel$summary.fixed), c("PA_intercept", "PO_intercept"))
  expect_equal(deparse1(RichModel$componentsJoint),
               "~-1 + shared_spatial(main = geometry, model = shared_field) + speciesShared(main = geometry, model = speciesField, group = speciesSpatialGroup, control.group = list(model = \"iid\", hyper = list(prec = list(prior = \"loggamma\", param = c(1, 5e-05))))) + PA_intercept(1) + PO_intercept(1) + speciesName_intercepts(main = speciesName, model = \"iid\", constr = TRUE, hyper = list(prec = list(fixed = TRUE, initial = log(INLA::inla.set.control.fixed.default()$prec))))")
  rm(RichModel)

})
