testthat::test_that('obtainCovariate can correctly obtain the covariate layer, and transform it to the desired projection', {

 skip_on_cran()

  #args
  covname <- "tavg"
  countries <- c('Norway', 'Sweden')
  projection <- '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
  #path <- './tests/testthat'
  path <- 'cov_file'
  cov <- obtainCovariate(covname, countries,
                         projection, path)

  expect_equal(class(cov)[1], 'SpatRaster')
  expect_identical(st_crs(cov)[2], st_crs(projection)[2])

  #Change CRS
  projection2 <- 'EPSG:4326'
  cov2 <- obtainCovariate(covname, countries,
                         projection2, path)

  expect_equal(class(cov2)[1], 'SpatRaster')
  expect_identical(st_crs(cov2)[2], st_crs(projection2)[2])

  unlink('cov_file', recursive = TRUE)



})
