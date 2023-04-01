obtainGBIF <- function(query,
                       datasetName,
                       gbifopts,
                       geometry,
                       projection,
                       varsKeep,
                       ...) {

  if (!inherits(geometry, 'sf')) geometry <- as(geometry, 'sf')

  boundaryCheck <- sf::st_transform(geometry, crs = CRS("+proj=longlat +ellps=WGS84"))

  speciesList <- vector(mode = 'list', length = length(query))
  names(speciesList) <- query

  for (species in query) {

  message(paste('Finding GBIF observations for:', species,'\n'))

  speciesOCC <- spocc::occ(query = species,
                        geometry = boundaryCheck,
                        gbifopts = gbifopts,
                        ...)

   if (nrow(data.frame(speciesOCC$gbif$data)) == 0) stop ('Species provided not available in specified area.')

  speciesSF <- st_as_sf(x = data.frame(speciesOCC$gbif$data[[1]]),
                        coords = c('longitude', 'latitude'),
                        crs = projection)

  speciesSF <- st_transform(speciesSF, crs = projection)

  speciesIn <- speciesSF[unlist(st_intersects(geometry, speciesSF)),]

  if (nrow(speciesIn) == 0) warning(paste(species, 'provided no occurrence reccords over the specified region.'))

  speciesList[[species]] <- speciesIn

  }

  speciesList <- do.call(rbind, speciesList); speciesList

  }
