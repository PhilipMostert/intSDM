testthat::test_that('sdmWorkflow produces the correct output given different Workflow situations.', {

  skip_on_cran()

  ##Create different workflows here:
   #1. Just GBIF data
  proj <- '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
  species <- c('Fraxinus excelsior', 'Arnica montana')
  workflow <- startWorkflow(Species = species,
                            saveOptions = list(projectName = 'testthatexample', projectDirectory = './tests'),
                            Projection = proj, Countries = 'Norway',
                            Quiet = TRUE, Save = TRUE)

  workflow$addGBIF(datasetName = 'GBIF_data') #Get less species
  expect_error(sdmWorkflow(Workflow = workflow)) #Test no output given
  workflow$workflowOutput('Model')
  expect_error(sdmWorkflow(Workflow = workflow)) #Test no mesh provided
  workflow$addMesh(max.edge = 200000)
  #Test something about CV-method -- none specified but given as output
  #Need to test a lot of the copy model; points spatial; points intercept parts

  workflow$addCovariates(worldClim = 'tmax')#add covariate here
  expect_true(dir.exists('./tests/testthatexample/Covariates'))

  sdmWorkflow(Workflow = workflow)
  expect_true(all(c(dir.exists('./tests/testthatexample/Fraxinus_excelsior'),
                    dir.exists('./tests/testthatexample/Arnica_montana'))))

  expect_true(all(c(file.exists('./tests/testthatexample/Fraxinus_excelsior/intModel.rds'),
                    file.exists('./tests/testthatexample/Arnica_montana/intModel.rds'))))

  Fraxinus_excelsior <- readRDS(file = './tests/testthatexample/Fraxinus_excelsior/intModel.rds')
  Fraxinus_excelsior
  #Test no output provided
  #Need to test that the model also saves all the outputs in a given directory, for multiple species.



   #2. GBIF + PA + Counts data
   #3. Change fields settings + other INLA options
   #4. Bias fields
   #5. Test all the different possible outcomes

})
