

formatStructured <- function(dataOCC, type, varsOld, varsNew,
                             coordinates, projection, boundary) {

  varsOld <- unlist(varsOld)
  varsNew <- unlist(varsNew)

  if (inherits(data, 'Spatial')) dataOCC <- as(dataOCC, 'sf')

  if (inherits(data, 'data.frame')) {

    if (nrow (dataOCC) <= 2) dataOCC$.__PORESP.__ <- 1

    dataOCC <- st_as_sf(x = dataOCC[, c('x', 'y', names(dataOCC)[names(dataOCC) %in% varsOld])],
                    coords = coordinates,
                    crs = projection)

  }

  initRows <- nrow(dataOCC)

  namesData <- colnames(dataOCC)[!colnames(dataOCC) %in% c('geometry', '.__PORESP.__')]

  oldNames <- varsOld[match(namesData, varsOld)]
  newNames <- varsNew[match(names(oldNames), names(varsNew))]

  colnames(dataOCC)[!colnames(dataOCC) %in% c('geometry', '.__PORESP.__')] <- newNames

  dataOCC <- dataOCC[unlist(st_intersects(boundary, dataOCC)),]

  if (nrow(dataOCC) == 0) stop('Dataset provided has no reccords over the boundary.')
  if (initRows < nrow(dataOCC)) warning('Some of the records provided are not over the boundary, and will therefore be removed.')

  dataOCC

}
