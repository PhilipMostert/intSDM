testthat::test_that('generateAbsences correctly creates absences for the data.', {

  ##First set up workflow
  skip_on_cran()
  skip_if_not(local_testthat_geodata_path())

  test_data <- readRDS(system.file('extdata/test_data.rds', package = 'intSDM'))

  POpoints <- list(Fraxinus_excelsior = test_data$POpoints)
  PApoints <- list(Fraxinus_excelsior = test_data$PApoints)
  Mesh <- test_data$Mesh

  proj <- '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
  species <- c('Fraxinus excelsior', 'Ulmus glabra', 'Arnica montana')

  workflow <- try(startWorkflow(Species = species,
                            saveOptions = list(projectName = 'testthatexample'),
                            Projection = proj,
                            Countries = c('Sweden', 'Norway'),
                            Quiet = TRUE, Save = FALSE))

  if (inherits(workflow, 'try-error')) {


    workflow <- startWorkflow(Species = species,
                  saveOptions = list(projectName = 'testthatexample'),
                  Projection = proj,
                  Quiet = TRUE, Save = FALSE)

    countries <- st_as_sf(geodata::world(path = geodata::geodata_path()))
    countries <- countries[countries$NAME_0 %in% c('Norway', 'Sweden'),]
    countries <- st_transform(countries, proj)

    workflow$addArea(Object = countries)

  }


  if (is.null(workflow$.__enclos_env__$private$Area)) {

    map <- st_as_sf(geodata::world(path = geodata::geodata_path()))
    map <- map[map$NAME_0 == 'Norway',]
    map <- st_transform(map, proj)

    workflow$addArea(Object = map)

  }

  for (spec in species[2:3]) {

    specName <- gsub(' ', '_', spec)

    POjitter <- st_jitter(POpoints$Fraxinus_excelsior, amount = 50)
    POjitter$name <- spec

    PAjitter <- st_jitter(PApoints$Fraxinus_excelsior, amount = 50)
    PAjitter$name <- spec

    POpoints[[specName]] <- POjitter
    PApoints[[specName]] <- PAjitter



  }

  POpointsComb <- do.call(rbind, POpoints)
  PApointsComb <- do.call(rbind, PApoints)

  workflow$addStructured(dataStructured = POpointsComb, datasetType = 'PO', datasetName = 'PO', speciesName = 'name')
  workflow$addStructured(dataStructured = PApointsComb, datasetType = 'PA', datasetName = 'PA', speciesName = 'name', responseName = 'pres', generateAbsences = FALSE)

  paData <- lapply(workflow$.__enclos_env__$private$dataStructured, function(x) x[['PA']])

  workflow$addStructured(dataStructured = PApointsComb, datasetType = 'PA', datasetName = 'PA', speciesName = 'name', responseName = 'pres', generateAbsences = TRUE)

  expect_true(all(names(workflow$.__enclos_env__$private$dataGBIF) %in% sub(" ", '_', species)))

  expect_true(all(unlist(lapply(workflow$.__enclos_env__$private$dataStructured, function(x) names(x))) %in% c('PO', 'PA')))

  expect_true(nrow(paData$Fraxinus_excelsior) < nrow(workflow$.__enclos_env__$private$dataStructured$Fraxinus_excelsior$PA))
  expect_true(nrow(paData$Ulmus_glabra) < nrow(workflow$.__enclos_env__$private$dataStructured$Ulmus_glabra$PA))
  expect_true(nrow(paData$Arnica_montana) < nrow(workflow$.__enclos_env__$private$dataStructured$Arnica_montana$PA))

  ##Test Richness = TRUE
  proj <- '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=km +no_defs'
  species <- c('Fraxinus excelsior', 'Ulmus glabra', 'Arnica montana')
  workflow <- try(startWorkflow(Species = species,
                                saveOptions = list(projectName = 'testthatexample'),
                                Projection = proj, Richness = TRUE,
                                Countries = c('Sweden', 'Norway'),
                                Quiet = TRUE, Save = FALSE))

  if (inherits(workflow, 'try-error')) {


    workflow <- startWorkflow(Species = species,
                              saveOptions = list(projectName = 'testthatexample'),
                              Projection = proj, Richness = TRUE,
                              Quiet = TRUE, Save = FALSE)

    countries <- st_as_sf(geodata::world(path = geodata::geodata_path()))
    countries <- countries[countries$NAME_0 %in% c('Norway', 'Sweden'),]
    countries <- st_transform(countries, proj)

    workflow$addArea(Object = countries)

  }



  if (is.null(workflow$.__enclos_env__$private$Area)) {

    map <- st_as_sf(geodata::world(path = geodata::geodata_path()))
    map <- map[map$NAME_0 == 'Norway',]
    map <- st_transform(map, proj)

    workflow$addArea(Object = map)

  }

  workflow$addStructured(dataStructured = POpointsComb, datasetType = 'PO', datasetName = 'PO', speciesName = 'name')
  workflow$addStructured(dataStructured = PApointsComb, datasetType = 'PA', datasetName = 'PA', speciesName = 'name', responseName = 'pres', generateAbsences = FALSE)

  paData <- lapply(workflow$.__enclos_env__$private$dataStructured, function(x) x[['PA']])

  #expect_warning(
   #{
      #workflow$addGBIF(datasetType = 'PA', datasetName = 'PA', generateAbsences = TRUE)
    #},
    #"datasetName already provided before. The older dataset will therefore be removed."
  #)
  workflow$addStructured(dataStructured = PApointsComb, datasetType = 'PA', datasetName = 'PA', speciesName = 'name', responseName = 'pres', generateAbsences = TRUE)

  expect_true(all(names(workflow$.__enclos_env__$private$dataStructured) %in% c('PO', 'PA')))

  expect_true(all(unlist(lapply(workflow$.__enclos_env__$private$dataGBIF, function(x) names(x))) %in% c('PO', 'PA')))

  expect_true(nrow(paData$PA[paData$PA$speciesName == 'Fraxinus_excelsior',]) < nrow(workflow$.__enclos_env__$private$dataStructured$PA$PA[workflow$.__enclos_env__$private$dataStructured$PA$PA$speciesName == 'Fraxinus_excelsior',]))
  expect_true(nrow(paData$PA[paData$PA$speciesName == 'Ulmus_glabra',]) < nrow(workflow$.__enclos_env__$private$dataStructured$PA$PA[workflow$.__enclos_env__$private$dataStructured$PA$PA$speciesName == 'Ulmus_glabra',]))
  expect_true(nrow(paData$PA[paData$PA$speciesName == 'Arnica_montana',]) < nrow(workflow$.__enclos_env__$private$dataStructured$PA$PA[workflow$.__enclos_env__$private$dataStructured$PA$PA$speciesName == 'Arnica_montana',]))



})
