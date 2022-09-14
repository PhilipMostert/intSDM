library(dplyr)
library(sp)
library(maps)
library(maptools)
library(RANN)

download_url_trd <- "https://gbif.vm.ntnu.no/ipt/archive.do?r=vascularplantfieldnotes"

if(!file.exists("Data/gbif_PA_trd.zip")) {
  download.file(url=download_url_trd, destfile="Data/gbif_PA_trd.zip", mode = "wb")
}
# Download data from University of Oslo
download_url_osl <- "https://ipt.gbif.no/archive.do?r=o_vxl"
if(!file.exists("Data/gbif_PA_osl.zip")) {
  download.file(url=download_url_osl, destfile="Data/gbif_PA_osl.zip", mode = "wb")
}
## 2.2. Read into R and prepare for filtering ----------
# Read NTNU data and select relevant columns
gbif_data_trd <- read.table(unzip(zipfile = "Data/gbif_PA_trd.zip", files="occurrence.txt"), header=T, sep="\t", quote="", fill=FALSE) %>%
  select(c(eventID, decimalLatitude, decimalLongitude,
           coordinateUncertaintyInMeters, scientificName,
           individualCount))
# Read data from UiO
gbif_data_oslo <- read.table(unzip(zipfile = "Data/gbif_PA_osl.zip", files="occurrence.txt"), header=T, sep="\t", quote="", fill=TRUE) %>%
  # remove missing coordinates
  subset(., !is.na(decimalLatitude) & !is.na(decimalLongitude)) %>%
  # remove data with missing information on year
  subset(., year != 0) %>%
  # make unique id for each event
  mutate(eventID = paste(day, month, year, decimalLatitude, decimalLongitude, sep = "_")) %>%
  # make column with information about presence for each event
  mutate(count = 1) %>%
  group_by(eventID, scientificName) %>%
  mutate(count = sum(count)) %>%
  mutate(individualCount = ifelse(count >= 1, 1, 0)) %>%
  # select relevant columns
  select(c(decimalLatitude, decimalLongitude,
           individualCount, scientificName, eventID, coordinateUncertaintyInMeters)) %>%
  data.frame()
## 2.3. Filter the PA data from the Red List Data species list (made in step 1) -----
## Will result in three files. (want visit numbers too, column per variable).
# download and read redlist data
download_url_redlist <- "https://artsdatabanken.no/Rodliste2015/sok/Eksport?kategori=re%2ccr%2cen%2cvu%2cnt%2cdd&vurderings%u00e5r=2015&vurderingscontext=n&taxonrank=species"
if(!file.exists("Data/Rodlista2015_Artsdatabanken_format.csv")) {
  download.file(url=download_url_redlist, destfile="Data/Rodlista2015_Artsdatabanken_format.csv", mode = "wb")
}
redlist <- read.csv2("Data/Rodlista2015_Artsdatabanken_format.csv", na.strings = c("", NA), fileEncoding="UTF-16LE")[c("Vitenskapelig.navn", "Kategori")]
# select species in PA files also present in redlist
# note: species names in PA file does not seem to be formatted in a consequent way,
# names formatted differently than the redlist will not be filtered
# one data frame for each relevant redlist category
#CR.PA.data.trd <- gbif_data_trd %>% subset(scientificName %in% redlist$Vitenskapelig.navn[redlist$Kategori %in% c("CR")])
#EN.PA.data.trd <- gbif_data_trd %>% subset(scientificName %in% redlist$Vitenskapelig.navn[redlist$Kategori %in% c("EN", "ENº")])
VU.PA.data.trd <- gbif_data_trd %>% subset(scientificName %in% redlist$Vitenskapelig.navn[redlist$Kategori %in% c("VU", paste0("VU", "\u00B0"))])
#CR.PA.data.osl <- gbif_data_oslo %>% subset(scientificName %in% redlist$Vitenskapelig.navn[redlist$Kategori %in% c("CR")])
#EN.PA.data.osl <- gbif_data_oslo %>% subset(scientificName %in% redlist$Vitenskapelig.navn[redlist$Kategori %in% c("EN", "ENº")])
#CR.PA.data <- rbind(CR.PA.data.trd, CR.PA.data.osl)
#EN.PA.data <- rbind(EN.PA.data.trd, EN.PA.data.osl)
VU.PA.data <- rbind(VU.PA.data.trd, VU.PA.data.osl)

##Make data spatial
Projection <- CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
norwayfill <- map("world", "norway", fill=TRUE, plot=FALSE,
                  ylim=c(58,72), xlim=c(4,32))
IDs <- sapply(strsplit(norwayfill$names, ":"), function(x) x[1])
norway.poly <- map2SpatialPolygons(norwayfill, IDs = IDs,
                                   proj4string = Projection)

abundant <- VU.PA.data %>% group_by(scientificName) %>%
  count() %>%
  arrange(desc(n)) %>%
  data.frame()
abundant <- abundant[1:3, 'scientificName']

PA_data <- VU.PA.data[VU.PA.data$scientificName%in%abundant,]

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
names(PA_redlist) <-c('species', 'individualCount')

usethis::use_data(PA_redlist, overwrite = TRUE, compress = 'xz')
