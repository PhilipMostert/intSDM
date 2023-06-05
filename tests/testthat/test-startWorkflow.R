testthat::test_that('startWorkflow can correctly create a species_model object, and store all the relavant metadata required by the workflow.', {

  skip_on_cran()

  expect_error(startWorkflow(), 'Please provide projectName in the saveOptions list.')
  expect_error(startWorkflow(saveOptions = list(projectName = 'testthat')), 'At least one species name needs to be provided.')

  countries <- c('Sweden', 'Norway')
  proj <- '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
  species <- 'Fraxinus excelsior'

  expect_error(startWorkflow(Countries = countries,
                            Species = species,
                            saveOptions = list(projectName = 'testthat')), 'argument "Projection" is missing, with no default')


  expect_message(startWorkflow(Countries = countries,
                Species = species,
                saveOptions = list(projectName = 'testthat'),
                Projection = proj,
                Quiet = TRUE), regexp = NA)

  expect_message(startWorkflow(Species = species,
                               saveOptions = list(projectName = 'testthat'),
                               Projection = proj,
                               Quiet = FALSE), regexp = NULL)

  workflow <- startWorkflow(Countries = countries,
                            Species = species,
                            saveOptions = list(projectName = 'testthatexample',
                                               projectDirectory = './'),
                            Projection = proj,
                            Quiet = TRUE)

  expect_setequal(class(workflow), c("species_model", "R6"))
  expect_true(dir.exists('./testthatexample'))

  unlink('./testthatexample', recursive = TRUE)


  })
