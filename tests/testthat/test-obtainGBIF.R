testthat::test_that('obtainGBIG can correctly obtain observations of species in the correct boundary, and transform it to the desired projection', {

  skip_on_cran()

  #args
  speciesIn <- 'Fraxinus excelsior'

  proj <- '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

  map <- obtainArea(names = c('Norway'), projection = proj)

  expect_message(obtainGBIF(query = speciesIn, country = 'NO',
                            datasettype = 'PO', geometry = map, projection = proj), 'Finding GBIF observations for: Fraxinus excelsior')

  species <- obtainGBIF(query = speciesIn, datasetName = 'dataset',
                        datasettype = 'PO', country = 'NO',
                        coordinateUncertaintyInMeters = 50,
                        geometry = map, projection = proj)

  expect_equal(class(species), c('sf', 'data.frame'))
  expect_true(all(species$coordinateUncertaintyInMeters <= 50))
  expect_identical(st_crs(species)[2], st_crs(proj)[2])

  speciesPA <- obtainGBIF(query = speciesIn, datasetName = 'dataset',
                          datasettype = 'PA',country = 'NO',
                          geometry = map, projection = proj)

  expect_setequal(unique(speciesPA$occurrenceStatus), c(0,1))

  speciesIn2 <- 'Ceratotherium simum'

  expect_error(obtainGBIF(query = speciesIn2, datasettype = 'PA',
                          country = 'NO',
                          geometry = map, projection = proj), 'Species provided not available in specified area.')

  ##Multiple years

  speciesTime <- obtainGBIF(query = speciesIn, datasetName = 'dataset',
                            datasettype = 'PO',country = 'NO',
                            geometry = map, projection = proj, year = 2010:2012)

  expect_setequal(unique(speciesTime$year), c(2010, 2011, 2012))

})