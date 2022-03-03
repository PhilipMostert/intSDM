testthat::test_that('Test that structured_data is able to differentiate between PO and PA data + provide correct output.', {

  #Arbitrary projection
  projection <- CRS('+proj=tmerc')
  #Make random shape to generate points on
  x <- c(16.48438,  17.49512,  24.74609, 22.59277, 16.48438)
  y <- c(59.736328125, 55.1220703125, 55.0341796875, 61.142578125, 59.736328125)
  xy <- cbind(x, y)

  Poly = Polygon(xy)
  Poly = Polygons(list(Poly),1)
  SpatialPoly = SpatialPolygons(list(Poly), proj4string = projection)

  #Make random points
  #Random presence only dataset
  PO <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  #Add species name
  PO$species <- sample(x = c('a','b'), size = nrow(PO@coords), replace = TRUE)

  #Random presence absence dataset
  PA <- spsample(SpatialPoly, n = 100, 'random', CRSobs = projection)
  #Add species name
  PA$species <- sample(x = c('c','d'), size = nrow(PA@coords), replace = TRUE)
  #Add response name
  PA$PAresp <- sample(x = c(0,1), size = nrow(PA@coords), replace = TRUE)
  #Add trial name
  PA$trials <- sample(x = c(1,2,3), size = nrow(PA@coords), replace = TRUE)

  dataObj <- structured_data(PO, PA, datasetType = c('PO','PA'), responsePA = 'PAresp',
                            trialsPA = 'trials', speciesName = 'species', coordinateNames = c('x','y'))

  expect_s4_class(dataObj, 'structuredData')

  expect_true(class(dataObj@dataPO$PO) == 'SpatialPointsDataFrame')
  expect_true(class(dataObj@dataPA$PA) == 'SpatialPointsDataFrame')

  expect_identical(attributes(dataObj)$datasetType, c('PO','PA'))
  expect_identical(attributes(dataObj)$responsePA, 'PAresp')
  expect_identical(attributes(dataObj)$trialsPA, 'trials')
  expect_identical(attributes(dataObj)$speciesName, 'species')
  expect_identical(attributes(dataObj)$coordinateNames, c('x','y'))

  ##Check warnings

  expect_error(structured_data(PO, PA, datasetType = c('Present only','PA'), responsePA = 'PAresp',
                                trialsPA = 'trials', speciesName = 'species', coordinateNames = c('x','y')),'datasetType must be a vector with values: "count", PO" or "PA".')
  PA_wrong_coords <- PA
  colnames(PA_wrong_coords@coords) <- c('NotX','NotY')
  expect_error(structured_data(PO, PA_wrong_coords, datasetType = c('PO','PA'), responsePA = 'PAresp',
                              trialsPA = 'trials', speciesName = 'species', coordinateNames = c('x','y')),'All datasets are required to have the same coordinate names specified by coordinateNames.')

  PO_wrong_species <- PO
  names(PO_wrong_species@data) <- 'NotSpecies'
  expect_error(structured_data(PO_wrong_species, PA, datasetType = c('PO','PA'), responsePA = 'PAresp',
                              trialsPA = 'trials', speciesName = 'species', coordinateNames = c('x','y')), 'All datasets names are required to have the same species variable name specified with speciesName.')


})
