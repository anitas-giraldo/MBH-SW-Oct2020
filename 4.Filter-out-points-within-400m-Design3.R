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



## read cluster centres ----
c1 <- readOGR(paste(d.dir, "InNOutMP-clusters-d6.shp", sep='/'))

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


# join with dean's actual points ----
cdf <- as.data.frame(c1)
str(cdf)
levels(cdf$class)
levels(cdf$class)[levels(cdf$class)=="Legacy"] <- "Deans"
cdf <- cdf[cdf$class!= "MBH",]
cdf <- droplevels(cdf)
coordinates(cdf) <- ~coords.x1+coords.x2

head(swbruvs)
head(cdf)
cdf <- cdf[,c(3,4,5)]
head(cdf)
proj4string(cdf) <- proj4string(swbruvs)


allb <- union(swbruvs, cdf)

# save --
writeOGR( allb, dsn=paste(d.dir, 'final', sep='/'), layer=paste( "InNOutMP-Bruvs-d3_utm", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR( allb, dsn=paste(d.dir, 'final', sep='/'), layer=paste( "InNOutMP-Bruvs-d3_latlong", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)






