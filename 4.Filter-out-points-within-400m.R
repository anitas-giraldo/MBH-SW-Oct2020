#### Get the clustered points design and filterout the points that are within 250 m of eachother ####

library( MBHdesign)
library( parallel)
library( class)
library( fields)
library( pdist)
library( raster)
library(rgeos)
library( rgdal)
library( sp)


# clear environment ----
rm(list = ls())


# Set working directory ####
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
p.dir <- paste(w.dir, "plots", sep ='/')
s.dir <- paste(w.dir, "spatial_data", sep='/')
#d.dir <- paste(w.dir, "Fish-HWY", sep='/')
d.dir <- paste(w.dir, "InNOutMP", sep='/')

## load raster of incl probs ----
ip <- raster(paste(d.dir, "inclProbs_InNOutMP-d3.tif", sep='/'))
plot(ip)

## read the spatial points created with MBH ----

p1 <- readOGR(paste(d.dir, "InNOutMP-BRUVS-d6.shp", sep='/')) # 240 features

proj4string(p1) <- proj4string(ip)

p2 <- readOGR(paste(d.dir, "FHWY-BRUVS-d3.shp", sep='/')) # 240 features

p3 <- readOGR("~/MBHdesignSW/design4/Bruvs4.3_2020-05-27.shp")

p4 <- readOGR("~/MBHdesignSW/design4/Bruvs4.4_2020-05-27.shp")

p5 <- readOGR("~/MBHdesignSW/design4/Bruvs4.5_2020-05-27.shp")

p6 <- readOGR("~/MBHdesignSW/design4/Bruvs4.6_2020-05-27.shp")

p8 <- readOGR("~/MBHdesignSW/design4/Bruvs4.8_2020-05-27.shp")

p9 <- readOGR("~/MBHdesignSW/design4/Bruvs4.9_2020-05-27.shp")

## read cluster centres ----
c1 <- readOGR("~/MBHdesignSW/design4/clusters4.1_2020-05-27.shp")

c2 <- readOGR("~/MBHdesignSW/design4/clusters4.2_2020-05-27.shp")

c3 <- readOGR("~/MBHdesignSW/design4/clusters4.3_2020-05-27.shp")

c4 <- readOGR("~/MBHdesignSW/design4/clusters4.4_2020-05-27.shp")

c5 <- readOGR("~/MBHdesignSW/design4/clusters4.5_2020-05-27.shp")

c6 <- readOGR("~/MBHdesignSW/design4/clusters4.6_2020-05-27.shp")

c7 <- readOGR("~/MBHdesignSW/design4/clusters4.7_2020-05-27.shp")

c8 <- readOGR("~/MBHdesignSW/design4/clusters4.8_2020-05-27.shp")

c9 <- readOGR("~/MBHdesignSW/design4/clusters4.9_2020-05-27.shp")


## read object in utm ----
ga <- raster(paste(s.dir, "ga4858_grid2_MSL.tiff", sep='/'))
gacrs <-  proj4string(ga)

## transform the points into UTM --
p1u <- spTransform(p1, gacrs)

## calculate if 2 points fall within 250m of eachother ----
# https://gis.stackexchange.com/questions/102796/remove-points-within-x-distance

## p1 ----
p1_matrix <- gWithinDistance(p1u, dist = 400, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p1_matrix[lower.tri(p1_matrix, diag=TRUE)] <- NA
p1_matrix

colSums(p1_matrix, na.rm=TRUE) == 0
v1 <- colSums(p1_matrix, na.rm=TRUE) == 0
p1[v1, ] # 98 features left

# plot --
plot(ip)
plot(p1[v1, ], pch=20, col="black", add=T)
#plot(c1, pch=20, col="blue", add=T)


# save points in utm and latlong ----

## Save points of design ---
swbruvs <- p1[v1, ]

# transform again
swbruvsu <- spTransform(swbruvs, gacrs)


# save --
writeOGR( swbruvsu, dsn=paste(d.dir, 'final', sep='/'), layer=paste( "InNOutMP-Bruvs-d3_utm", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR( swbruvs, dsn=paste(d.dir, 'final', sep='/'), layer=paste( "InNOutMP-Bruvs-d3_latlong", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)



## p2 ----
p2_matrix <- gWithinDistance(p2, dist = 250, byid = TRUE)
diag(p2_matrix) <- NA
p2_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p2_matrix[lower.tri(p2_matrix, diag=TRUE)] <- NA
p2_matrix

colSums(p2_matrix, na.rm=TRUE) == 0
v2 <- colSums(p2_matrix, na.rm=TRUE) == 0
p2[v2, ] # 14 features left

# plot --
plot(ip)
plot(p2[v2, ], pch=20, col="black", add=T)
plot(c2, pch=20, col="blue", add=T)

## p3 ----
p3_matrix <- gWithinDistance(p3, dist = 250, byid = TRUE)
diag(p3_matrix) <- NA
p3_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p3_matrix[lower.tri(p3_matrix, diag=TRUE)] <- NA
p3_matrix

colSums(p3_matrix, na.rm=TRUE) == 0
v3 <- colSums(p3_matrix, na.rm=TRUE) == 0
p3[v3, ] # 15 features left

# plot --
plot(ip, main="Opt1 - Clusters = 4, Points = 15")
plot(p3[v3, ], pch=20, col="black", add=T)
plot(c3, pch=20, col="blue", add=T)



## p4 ----
p4_matrix <- gWithinDistance(p4, dist = 250, byid = TRUE)
diag(p4_matrix) <- NA
p4_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p4_matrix[lower.tri(p4_matrix, diag=TRUE)] <- NA
p4_matrix

colSums(p4_matrix, na.rm=TRUE) == 0
v4 <- colSums(p4_matrix, na.rm=TRUE) == 0
p4[v4, ] # 12 features left

# plot --
plot(ip)
plot(p4[v4, ], pch=20, col="black", add=T)
plot(c4, pch=20, col="blue", add=T)


## p5 ----
p5_matrix <- gWithinDistance(p5, dist = 250, byid = TRUE)
diag(p5_matrix) <- NA
p5_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p5_matrix[lower.tri(p5_matrix, diag=TRUE)] <- NA
p5_matrix

colSums(p5_matrix, na.rm=TRUE) == 0
v5 <- colSums(p5_matrix, na.rm=TRUE) == 0
p5[v5, ] # 15 features left

# plot --
plot(ip)
plot(p5[v5, ], pch=20, col="black", add=T)
plot(c5, pch=20, col="blue", add=T)


## p6 ----
p6_matrix <- gWithinDistance(p6, dist = 250, byid = TRUE)
diag(p6_matrix) <- NA
p6_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p6_matrix[lower.tri(p6_matrix, diag=TRUE)] <- NA
p6_matrix

colSums(p6_matrix, na.rm=TRUE) == 0
v6 <- colSums(p6_matrix, na.rm=TRUE) == 0
p6[v6, ] # 13 features left

# plot --
plot(ip)
plot(p6[v6, ], pch=20, col="black", add=T)
plot(c6, pch=20, col="blue", add=T)


## p7 ----
p7_matrix <- gWithinDistance(p7, dist = 250, byid = TRUE)
diag(p7_matrix) <- NA
p7_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p7_matrix[lower.tri(p7_matrix, diag=TRUE)] <- NA
p7_matrix

colSums(p7_matrix, na.rm=TRUE) == 0
v7 <- colSums(p7_matrix, na.rm=TRUE) == 0
p7[v7, ] # 12 features left

# plot --
plot(ip, main = " Opt2 - Clusters = 3, Points = 12")
plot(p7[v7, ], pch=20, col="black", add=T)
plot(c7, pch=20, col="blue", add=T)


## p8 ----
p8_matrix <- gWithinDistance(p8, dist = 250, byid = TRUE)
diag(p8_matrix) <- NA
p8_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p8_matrix[lower.tri(p8_matrix, diag=TRUE)] <- NA
p8_matrix

colSums(p8_matrix, na.rm=TRUE) == 0
v8 <- colSums(p8_matrix, na.rm=TRUE) == 0
p8[v8, ] # 10 features left

# plot --
plot(ip)
plot(p8[v8, ], pch=20, col="black", add=T)
plot(c8, pch=20, col="blue", add=T)


## p9 ----
p9_matrix <- gWithinDistance(p9, dist = 250, byid = TRUE)
diag(p9_matrix) <- NA
p9_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p9_matrix[lower.tri(p9_matrix, diag=TRUE)] <- NA
p9_matrix

colSums(p9_matrix, na.rm=TRUE) == 0
v9 <- colSums(p9_matrix, na.rm=TRUE) == 0
p9[v9, ] # 12 features left

# plot --
plot(ip)
plot(p9[v9, ], pch=20, col="black", add=T)
plot(c9, pch=20, col="blue", add=T)


#### Chosen design 7 ####

library(viridis)

pdf(paste(w.dir, "Fish-HWY", "FishHWY-Bathy-Map-design5.pdf", sep='/'), height=7, width=8)
plot(b, col=viridis(100), main = "Recreational fishermen highway sampling design")
plot(p1[v1, ], pch=20, cex=1, col="black", add=T)
plot(p1[v1, ], pch=20, cex=0.6, col="white", add=T)
#plot(c7, pch=20, col="blue", add=T)
dev.off()

## Save points of design 7 ----
swbruvs <- p7[v7, ]
swbruvs

## transform points to Lat long --
ras <- raster("~/MBHdesignGB/SpatialData/GB_CMR_bathy.tif")
ras

td <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

swbruvst <- spTransform(swbruvs, td)


writeOGR( swbruvs, dsn="~/MBHdesignSW/design4", layer=paste( "SW-Bruvs_utm", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR( swbruvst, dsn="~/MBHdesignSW/design4", layer=paste( "SW-Bruvs_latlong", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)





