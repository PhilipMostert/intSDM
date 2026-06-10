library(sf)
library(INLA)
library(terra)
library(giscoR)
library(inlabru)
library(fmesher)

set.seed(123)

crs <- '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=km +no_defs'

Countries <- giscoR::gisco_countries_2024

NorwaySweden <- Countries[Countries$NAME_ENGL %in% c('Norway', 'Sweden'),]

NorwaySweden <- st_transform(NorwaySweden, crs)

NorwaySweden <- st_buffer(NorwaySweden, 50)

rastNS2 <- rastNS1 <- rast(NorwaySweden, nrows = 100, ncols = 100)

y0 <- x0 <- seq(0, 10, length.out = ceiling(sqrt(nrow(terra::crds(rastNS1)))))

values(rastNS1) <- (outer(y0,x0,
                            function (x,y) 1/2*x - 2*y)/3)[1:nrow(values(rastNS1))]

values(rastNS2) <- (outer(y0,x0,
                          function (x,y) 1/2*x + 1/5*y - 5)/2)[1:nrow(values(rastNS2))]

rastNS <- c(rastNS1, rastNS2)

names(rastNS) <- c('Cov1', 'Cov2')

#Coarse mesh
Mesh <- fmesher::fm_mesh_2d(boundary = st_union(NorwaySweden),
                            max.edge = c(50, 150) * 2,
                            cutoff = 30,
                            offset = c(100, 500), crs = crs)


PApoints <- as.data.frame(rastNS$Cov1, xy = TRUE,
                          cells = FALSE)

POpoints <- as.data.frame(rastNS$Cov2, xy = TRUE,
                          cells = FALSE)

PApoints <- st_as_sf(x = PApoints,
               coords = c('x', 'y'), crs = crs)

POpoints <- st_as_sf(x = POpoints,
                     coords = c('x', 'y'), crs = crs)

PApoints <- st_filter(PApoints, NorwaySweden)

POpoints <- st_filter(POpoints, NorwaySweden)

POpoints$keep <- rbinom(n = nrow(POpoints), size =  1, prob = INLA::inla.link.logit(x = POpoints$Cov2, inverse = TRUE))
POpoints <- POpoints[POpoints$keep == 1,]

PApoints <- PApoints[sample.int(n = nrow(PApoints), size = 100 , replace = FALSE),]

PApoints$pres <- rbinom(n = nrow(PApoints), size = 1, prob = INLA:::inla.link.cloglog(x = 0 + PApoints$Cov1, inverse = TRUE))


POpoints$name <- PApoints$name <- 'Fraxinus excelsior'

test_data <- list(PApoints = PApoints,
                      POpoints = POpoints,
                      Mesh = Mesh)

terra::writeRaster(
  rastNS,
  filename = 'inst/extdata/test_covariates.tif',
  overwrite = TRUE,
  gdal = c("COMPRESS=LZW")
)

saveRDS(test_data, file = 'inst/extdata/test_data.rds')

#usethis::use_data(test_data, overwrite = TRUE)

