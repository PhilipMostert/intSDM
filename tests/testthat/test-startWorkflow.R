testthat::test_that('startWorkflow can correctly create a species_model object, and store all the relavant metadata required by the workflow.', {

  skip_on_cran()
  skip_if_not(local_testthat_geodata_path())

  projectDir <- tempfile("intSDM_startWorkflow_test_")
  on.exit(unlink(projectDir, recursive = TRUE))
  dir.create(projectDir, recursive = TRUE, showWarnings = FALSE)

  expect_error(startWorkflow(), 'Please provide projectName in the saveOptions list.')
  expect_error(startWorkflow(saveOptions = list(projectName = 'testthat',
                                                projectDirectory = projectDir)),
               'At least one species name needs to be provided.')

  countries <- c('Sweden', 'Norway')
  proj <- '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
  species <- 'Fraxinus excelsior'

  expect_error(startWorkflow(Countries = countries,
                             Species = species,
                             saveOptions = list(projectName = 'testthat',
                                                projectDirectory = projectDir)),
               'argument "Projection" is missing, with no default')

  countriesTry <- try(giscoR::gisco_countries_2024[giscoR::gisco_countries_2024$NAME_ENGL %in% c('Sweden', 'Norway'), ])

  skip_if(inherits(countriesTry, 'try-error'))

  expect_message(startWorkflow(Species = species,
                               saveOptions = list(projectName = 'testthat',
                                                  projectDirectory = projectDir),
                               Projection = proj,
                               Quiet = TRUE), regexp = NA)

  expect_message(startWorkflow(Species = species,
                               saveOptions = list(projectName = 'testthat',
                                                  projectDirectory = projectDir),
                               Projection = proj,
                               Quiet = FALSE), regexp = NULL)

  workflow <- startWorkflow(Species = species,
                            saveOptions = list(projectName = 'testthatexample',
                                               projectDirectory = projectDir),
                            Projection = proj,
                            Quiet = TRUE)

  expect_setequal(class(workflow), c("species_model", "R6"))
  expect_true(dir.exists(file.path(projectDir, 'testthatexample')))

  ##Test Richness model
  workflow <- startWorkflow(Species = species,
                            saveOptions = list(projectName = 'testthatexample',
                                               projectDirectory = projectDir),
                            Projection = proj, Richness = TRUE,
                            Save = FALSE,
                            Quiet = TRUE)
  expect_true(workflow$.__enclos_env__$private$richnessEstimate)


})
