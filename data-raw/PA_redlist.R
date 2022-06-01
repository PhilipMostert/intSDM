## code to prepare `DATASET` dataset goes here
library(sp)
library(maps)
library(maptools)
library(RANN)
Projection <- CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

norwayfill <- map("world", "norway", fill=TRUE, plot=FALSE,
                  ylim=c(58,72), xlim=c(4,32))
IDs <- sapply(strsplit(norwayfill$names, ":"), function(x) x[1])
norway.poly <- map2SpatialPolygons(norwayfill, IDs = IDs,
                                   proj4string = Projection)


norge_data <- read.csv('./data-raw/data/VU.PA.data.Norge.txt', sep = ' ')

olso_data <- read.csv('./data-raw/data/VU.PA.data.Olso.txt', sep = ' ')
olso_data <- olso_data[!(olso_data$decimalLongitude == 0),]
oslo_data <- olso_data[!(olso_data$decimalLatitude == 0),]

PA_data <- rbind(norge_data, olso_data)
abundant <- PA_data %>% group_by(scientificName) %>%
  count() %>%
  arrange(desc(n)) %>%
  data.frame()
abundant <- abundant[1:3, 'scientificName']

PA_data <- PA_data[PA_data$scientificName%in%abundant,]

PA_data <- sp::SpatialPointsDataFrame(coords = data.frame(PA_data$decimalLongitude, PA_data$decimalLatitude),
                                      data = data.frame(scientificName = PA_data$scientificName,
                                                        individualCount = PA_data$individualCount),
                                      proj4string = Projection)

grid <- makegrid(norway.poly, cellsize = 0.25, pretty = FALSE)
grid <- SpatialPoints(grid, proj4string = Projection)

spgrdWithin <- SpatialPixels(grid[norway.poly,])
spgrdWithin <- as(spgrdWithin, "SpatialPolygons")
PA_data_frame <- as.data.frame(PA_data)
closest <- RANN::nn2(data = grid@coords, query = PA_data_frame[,3:4], k = 1)
grid$ID <- seq(nrow(grid@coords))

closest_ind <- as.data.frame(closest) %>%
  dplyr::rename(ID = nn.idx)

closest_ind$scientificName <- PA_data_frame$scientificName
closest_ind$individualCount <- PA_data_frame$individualCount

joined <- dplyr::inner_join(closest_ind,as.data.frame(grid), by = 'ID')

grid_data <- sp::SpatialPointsDataFrame(coords = data.frame(joined$x1, joined$x2),
                                        data = data.frame(scientificName = joined$scientificName, individualCount = joined$individualCount),
                                        proj4string = Projection)
colnames(grid_data@coords) <- c('x','y')

#Remove the observations not in the boundary
grid_data <- grid_data[!is.na(over(grid_data, norway.poly)),]


species <- unique(grid_data$scientificName)

unique_grid <- as.data.frame(grid_data) %>%
  group_by(x,y) %>%
  slice(1) %>%
  dplyr::select(x,y) %>%
  data.frame()

grid_index <- list()
for (i in 1:nrow(unique_grid)) {

  x <- unique_grid[i,1]
  y <- unique_grid[i,2]

  present <- grid_data[grid_data@coords[,1] == x & grid_data@coords[,2] == y,]

  not_in <- species[!species%in%present$scientificName]

  if (!identical(not_in, character(0))) {

    ind <- data.frame(scientificName = not_in, individualCount = 0)

    ind_coords <- data.frame(x,y) %>% slice(rep(1:n(), each = length(not_in)))

    absent <- sp::SpatialPointsDataFrame(coords = ind_coords, data = ind,
                                         proj4string = Projection)

    grid_index[[i]] <- rbind.SpatialPointsDataFrame(present, absent)

  }
  else grid_index[[i]] <- present

}
PA_redlist <- do.call(rbind.SpatialPointsDataFrame, grid_index)
colnames(PA_redlist@coords) <- c('longitude','latitude')
PA_redlist$scientificName <- gsub(' ','_',PA_redlist$scientificName)
names(PA_redlist) <- c('species', 'individualCount')

usethis::use_data(PA_redlist, overwrite = TRUE, compress = 'xz')
