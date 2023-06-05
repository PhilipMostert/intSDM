testthat::test_that('obtainCovariate can correctly obtain the covariate layer, and transform it to the desired projection', {

 skip_on_cran()

  #args
  covname <- "tavg"
  countries <- c('Norway', 'Sweden')
  projection <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
  #path <- './tests/testthat'
  path <- './wc2.1_country'
  dir.create(path)

  if(dir.exists(path)) {

    ##If taking too long//servers down

  try(cov <- R.utils::withTimeout(
    obtainCovariate(covname, countries,
                    projection, path), timeout = 60, onTimeout = 'silent'))

  if ('cov' %in% ls()) {

  expect_equal(class(cov)[1], 'SpatRaster')
  expect_identical(st_crs(cov)[2], st_crs(projection)[2])


  }
  #Change CRS
  projection2 <- 'EPSG:4326'

  try(cov2 <-  R.utils::withTimeout(obtainCovariate(covname, countries,
                         projection2, path), timeout = 60, onTimeout = 'silent'))

  if ('cov2' %in% ls()) {

  expect_equal(class(cov2)[1], 'SpatRaster')
  expect_identical(st_crs(cov2)[2], st_crs(projection2)[2])
}
  unlink(path, recursive = TRUE)


  }

})
