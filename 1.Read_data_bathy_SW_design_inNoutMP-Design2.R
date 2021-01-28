### Read data ####

library( rgdal)
library( sp)
library( raster)


# Set working directory ####
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
s.dir <- paste(w.dir, "spatial_data", sep ='/')


# read poly of marine parks ----
mp <- readOGR(paste(s.dir, "SW_CMR_NP.shp", sep='/'))
mp$ORIG_NAME@data

# make polygon of area outside MP and inside ----

## Outside mp poly----
xcoord <- c(114.9312,114.9312,114.45724,114.45686,114.9312)
ycoord <- c(-34.0111,-34.0498,-34.0499,-34.01680,-34.0111)
xym <- cbind(xcoord, ycoord)
xym

pout <- Polygon(xym)
ps <- Polygons(list(pout),1)
poly.out <- SpatialPolygons(list(ps))
proj4string(poly.out) <- proj4string(mp)
# make spatial points data frame
data1 <- data.frame(f='Out')
poly.out.df <- SpatialPolygonsDataFrame(poly.out, data1)
plot(poly.out.df)


## Inside mp poly----
xcoord <- c(114.9314,114.9556,114.4553,114.45729,114.9314)
ycoord <- c(-34.0507,-34.1505,-34.1504,-34.05035,-34.0507)
xym <- cbind(xcoord, ycoord)
xym

pout <- Polygon(xym)
ps <- Polygons(list(pout),1)
poly.in <- SpatialPolygons(list(ps))
proj4string(poly.in) <- proj4string(mp)
# make spatial points data frame
data2 <- data.frame(f='In')
poly.in.df <- SpatialPolygonsDataFrame(poly.in, data2)
plot(poly.in.df)

#########################
#read in the country boundary -- for plotting mostly
#coastLine <- readOGR( dsn="C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/nsaasr9nnd_02211a04es_geo___ (1)")
#proj4string( coastLine) <- CRS("+init=epsg:4283")
#coastLine <- spTransform( coastLine, "+init=epsg:4326")

############################
#read in survey areas data
# zones is a list of polygons

writeOGR(poly.in.df, s.dir, "PolyInMP-d2", driver="ESRI Shapefile")
writeOGR(poly.out.df, s.dir, "PolyOutMP-d2", driver="ESRI Shapefile")

zones <- list()
zones$InsideMP <- poly.in.df
zones$OutsideMP <- poly.out.df
zones$SWMP <- mp
# join both polygons
zones$All <- raster::union(zones$SWMP, zones$OutsideMP)
zones$All2 <- raster::union(zones$InsideMP, zones$OutsideMP)
plot(zones$All2)



#intial look to see area
plot( zones$All, border='black')
plot( zones$InsideMP, add=TRUE, col='orange')
plot( zones$OutsideMP, add=TRUE, col='green')

# Read bathy ----
b <- raster(paste(s.dir, "GB-SW_250mBathy.tif", sep='/'))
plot(b)
plot(zones$All2, add=T)
# crop to area --
b2 <- crop(b, zones$All2)
plot(b2)
# remove cells deeper than 240m --
values(b2)[values(b2) < -270 ] = NA
# remove cells shallower than 35m --
values(b2)[values(b2) > -35 ] = NA

plot(b2)
plot(mp, add=T)

# save new bathy ----
writeRaster(b2, paste(s.dir, "SWBathy-InNOutMP-D2.tif", sep='/'), overwrite=T)

## read points from Dean ----
wp <- read.csv(paste(w.dir, "SW-points-Dean.csv", sep='/'))
# make sp --
coordinates(wp) <- ~Longitude+Latitude
points(wp)
proj4string(wp) <- proj4string(mp)

## calculate TPI ----
tpi <- terrain(b2, 'tpi')
plot(tpi)
points(wp)

writeRaster(tpi, paste(s.dir, "SW-TPI-InNOutMP-D2.tif", sep='/'), overwrite=T)

# calculate Slope ----
slope <- terrain(b2,'slope')
plot(slope)
points(wp)

writeRaster(slope, paste(s.dir, "SW-Slope-InNOutMP-D2.tif", sep='/'), overwrite =T)

# calculate aspect ----
aspect <- terrain(b2, 'aspect')
plot(aspect)
points(wp)




###################################################
#### converting polygons to a common raster.

#r <- gb_rasters$bathy

#plot( extent( r), add=TRUE)
#survArea first
###       ###       ### this takes a while for fine res data  ###      ###       ###
in_raster <- rasterize( x=zones$InsideMP, y=b2, field=zones$InsideMP@data[,1], bkg.value=-999, fun="first")
plot(in_raster)
out_raster <- rasterize( zones$OutsideMP, y=b2, field=zones$OutsideMP@data[,1], bkg.value=-999, fun="first")
plot(out_raster)

###################################
#convert and combine
tmp1 <- as.data.frame( in_raster, xy=TRUE)
tmp2 <- as.data.frame( out_raster, xy=TRUE)
tmp3 <- as.data.frame( b2, xy=TRUE)
tmp4 <- as.data.frame( slope, xy=TRUE)
tmp5 <- as.data.frame( tpi, xy=TRUE)


SWDat <- cbind( tmp1, tmp2[,3])
SWDat <- cbind( SWDat, tmp3[,3])
SWDat <- cbind( SWDat, tmp4[,3])
SWDat <- cbind( SWDat, tmp5[,3])
head(SWDat)

colnames( SWDat) <- c("Eastern", "Northing", "InsideMP", "OutsideMP", "Bathy", "Slope", "Tpi")


## New directory ----
d.dir <- paste(w.dir, "InNOutMP", sep='/')
setwd(d.dir)

saveRDS( SWDat, file="SWData_forInNOutMP-2.RDS")

sw_rasters <- list()
sw_rasters$bathy <- b2
sw_rasters$slope <- slope
sw_rasters$tpi <- tpi
saveRDS( sw_rasters, file="SWRasters_forInNOutMP-d2.RDS")
saveRDS(zones, file="SWZones_forInNOutMP-d2.RDS")

#rm( MUZ_raster, SPZ_raster, HPZ_raster, r, NPZ_raster, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)

#gc()


### for historical sites ----
plot(b2)
plot(zones$All2, add=T)
points(wp, add=T)

# Points in or out marine park ----
Dpoints <- raster::extract(zones$All2, wp, sp=T)
Dpoints

wpdf <- as.data.frame(wp)

Dpoints2 <- cbind(wpdf, Dpoints)
head(Dpoints2)
# remove the points outside --
Dpoints2<- Dpoints2[!is.na(Dpoints2$f),]
str(Dpoints2) #43 obs
coordinates(Dpoints2) <- ~Longitude+Latitude

plot(b2)
plot(zones$All2, add=T)
points(Dpoints2, col=Dpoints2$f, pch=20)

# Depth of points ----
Dpoints3 <- raster::extract(b2, Dpoints2, df=T)
head(Dpoints3)
str(Dpoints3)

df2 <- as.data.frame(Dpoints2)
Dpoints4 <- cbind(df2, Dpoints3)
head(Dpoints4) 
# remove points without bathy --
Dpoints4<- Dpoints4[!is.na(Dpoints4$GB.SW_250mBathy),]
str(Dpoints4) # 42 obs.
Dpoints5 <- Dpoints4
coordinates(Dpoints5) <- ~Longitude+Latitude

plot(b2)
plot(zones$All2, add=T)
points(Dpoints5, col=Dpoints5$f, pch=20)




write.csv(Dpoints5, paste(d.dir, "DeansPoints_forinNOutMP-d2.csv", sep ='/'))
BRUVS <- SpatialPointsDataFrame( coords=Dpoints5[,c("Longitude","Latitude")], data=Dpoints5, 
                                 proj4string = CRS( proj4string( mp)))
saveRDS( Dpoints5, file="DeansPoints_forinNOutMP-d2.RDS")




