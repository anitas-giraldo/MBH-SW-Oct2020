#install.packages("measurements")
library(measurements)

# Set working directory ####
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
s.dir <- paste(w.dir, "spatial_data", sep ='/')

# read data ----
d<-read.csv(paste(w.dir, "SW-coordinates_from_Dean.csv", sep='/'))
d
head(d)
str(d)

ds <-d

# add the negative sign to the Latitude
ds$Latitude <- paste('-', ds$Latitude, sep ='')
head(ds)

# convert from decimal minutes to decimal degrees
ds$Latitude = measurements::conv_unit(ds$Latitude, from = 'deg_dec_min', to = 'dec_deg')
ds$Longitude = measurements::conv_unit(ds$Longitude, from = 'deg_dec_min', to = 'dec_deg')

head(ds)

write.csv(ds, paste(w.dir, "SW-points-Dean.csv", sep='/'))


coordinates(ds) <- ~Latitude+Longitude

ds
bathy
proj4string(ds) <- proj4string(bathy)
ds2 <- spTransform(ds, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
ds2

writeOGR( ds2, dsn="~/MBHdesignGB/Design3", layer=paste( "hires_points", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)


head(dsp2)

# change the degree symbol to a space
dsp2$y = gsub('.', ' ', dsp2$y)
dsp2$x = gsub('.', ' ', dsp2$x)

# convert from decimal minutes to decimal degrees
dsp2$lat = measurements::conv_unit(dsp2$y, from = 'dec_deg', to = 'deg_dec_min')
dsp2$long = measurements::conv_unit(dsp2$x, from = 'dec_deg', to = 'deg_dec_min')



lr <- read.csv("~/MBHdesignGB/Design4_notClusTrans_LatLon (1).csv")
head(lr)
# convert from decimal minutes to decimal degrees
lr$y.1= measurements::conv_unit(lr$y.1, from = 'dec_deg', to = 'deg_dec_min')
lr$x.1 = measurements::conv_unit(lr$x.1, from = 'dec_deg', to = 'deg_dec_min')

write.csv(lr, "coordinates_lores.csv")

lrs <- lr

bathy
b

coordinates(lrs) <- ~x+y
proj4string(lrs) <- proj4string(bathy)
lrs2 <- spTransform(lrs, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
lrs2

writeOGR( lrs2, dsn="~/MBHdesignGB/Design3", layer=paste( "lowres_points", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)

ldf <- as.data.frame(lrs2)
head(ldf)
