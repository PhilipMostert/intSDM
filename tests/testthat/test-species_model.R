testthat::test_that('Test that species_model can produce all of the returns correctly.', {
  skip_on_cran()

  ##Create arbitrary data
  #Arbitrary projection
  projection <- CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

  ##Run species_model

  #Check return = 'boundary'

  boundary <- species_model(return = 'boundary', projection = projection)
  expect_true(class(boundary)[1] == 'SpatialPolygonsDataFrame')

  mesh <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(boundary),
                             cutoff=0.08, max.edge=c(1, 3), offset=c(1,1), crs = projection)


  #Check different boundary

  boundaryHedmark <- species_model(return = 'boundary', location = 'Hedmark')
  expect_true(class(boundaryHedmark)[1] == 'SpatialPolygonsDataFrame')
  expect_true(boundaryHedmark$NAME_1 == 'Hedmark')

  #Check an error

  expect_error(species_model(return = 'boundary', location = 'Sweden'),
               'At least one of the locations provided is not a valid county in Norway. NOTE: Trøndelag is given as Nord-Trøndelag and Sør-Trøndelag')

  ##Make species

  #Make random points
  #Random presence only dataset
  #Choose 10 points of arbitrary species
  species <- c('Fraxinus excelsior', 'Ulmus glabra')

  PO <- spsample(boundary, n = 20, 'random', CRSobs = projection)
  #Add species name
  PO$species <- sample(x = species, size = nrow(PO@coords), replace = TRUE)

  #Random presence absence dataset
  PA <- spsample(boundary, n = 20, 'random', CRSobs = projection)
  #Add species name
  PA$species <- sample(x = species, size = nrow(PA@coords), replace = TRUE)
  #Add response name
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA@coords), replace = TRUE)
  #Add trial name
  PA$trials <- sample(x = c(1,2,3), size = nrow(PA@coords), replace = TRUE)

  dataObj <- structured_data(PO, PA, datasetType = c('PO','PA'), responsePA = 'PAresp',
                             trialsPA = 'trials', speciesName = 'species', coordinateNames = c('x','y'))


  #Check return = 'species plot'

  speciesPlot <- species_model(return = 'species plot', structuredData = dataObj,
                               boundary = boundary, speciesNames = species,
                               limit = 10, mesh = mesh)

  expect_setequal(class(speciesPlot), c('gg', 'ggplot'))

  #Check return = 'model'
  if (requireNamespace('INLA')) {

  skip('Needs too much time to run.')

  covariate <- 'Annual Mean Temperature'
  model <- species_model(return = 'model',
                         boundary = boundary, speciesNames = species,
                         limit = 10, mesh = mesh, worldclimCovariates = covariate)

  expect_setequal(class(model), c('bruSdm', 'bru', 'iinla', 'inla'))

  #Check return = 'predictions'

  predictions <- species_model(return = 'predictions',
                               boundary = boundary, speciesNames = species,
                               limit = 10, mesh = mesh, worldclimCovariates = covariate)

  expect_setequal(class(predictions), c("predict_bruSdm", "list"))

  #Check return = 'predictions map'

  predictions <- species_model(return = 'predictions map',
                               boundary = boundary, speciesNames = species,
                               limit = 10, mesh = mesh, worldclimCovariates = covariate)

  expect_setequal(class(predictions), c('gg', 'ggplot'))
 }



})
