### Create inclusion probabilities ####

## for multibeam and lidar #######
#install.packages("MBHdesign")
library( MBHdesign)
library( parallel)
library( class)
library( fields)
#install.packages("pdist")
library( pdist)
library( raster)
library( rgdal)
library( sp)



# Set working directory ####
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
s.dir <- paste(w.dir, "spatial_data", sep ='/')


#########################
#read in data -----

#GBDat <- readRDS( "GBData_forDesign3.RDS")
#gb_rasters <- readRDS("GBRasters_forDesign3.RDS")
#zones <- readRDS( "GBZones_forDesign3.RDS")

b <- raster(paste(s.dir, "SWbathy_fishHWY.tif", sep='/'))
plot(b)
# remove all cells shallower than 35m
values(b)[values(b) > -35 ] = NA
plot(b)
# remove all cells deeper than 60m
values(b)[values(b) < -60 ] = NA
plot(b)

swnp <- readOGR(paste(s.dir, "SW_CMR_NP.shp", sep='/'))
plot(swnp, add=T)

## crop to desired location
#e <- drawExtent()
#b2 <- crop(b2, e)
#plot(b2)
#plot(swnp, add=T)
# save new raster
#writeRaster(b, paste(s.dir, "SWbathy_fishHWY.tif", sep ='/'), overwrite=T)

## read points from Dean ----
wp <- read.csv(paste(w.dir, "SW-points-Dean.csv", sep='/'))
# make sp --
coordinates(wp) <- ~Longitude+Latitude
points(wp)


####################################
####  Straw man for numbers of samples in each region
####################################

#straw.nums <- c( 16, 12, 6, 9)  # numbers of drops
#straw.props <- straw.nums / sum( straw.nums)
#names( straw.nums) <- names( straw.props) <- c( "MUZ", "SPZ", "HPZ", "NPZ")
#saveRDS( straw.nums, file="StrawmanNumbers_Zones.RDS")

#setwd("~/MBHdesignGB/Design3/")


###################################
####  Hand-picking Bathy cut points
####  And their numbers of drops
###################################

#Bathy.quant <- c(0,0.8,0.9,0.925,0.95,0.975,1)
Bathy.quant <- c(0,0.4,0.6,0.7,1)
Bathy.cuts <- quantile(b, Bathy.quant)#c( -Inf,0.02,0.04,0.08,0.16,Inf)
#trying to make it so there is no hand-picking (except for the hand-picked function)
tmp <- cumsum(Bathy.quant)
Bathy.targetNums <- rep(floor(18/8), 4)#floor( ( tmp / sum( tmp))[-1] * 200)#rep( 40, 5)#c( 20,20,30,65,65)
Bathy.targetProps <-  Bathy.targetNums / sum( Bathy.targetNums)

########################
#depth limiting and some elementary cleaning

#minDepthInAMP <- max( NingalooDat$BATHY, na.rm=TRUE)

#NingalooDat[ is.na( NingalooDat$BATHY) | NingalooDat$BATHY < -195 & NingalooDat$BATHY > minDepthInAMP, c("BATHY","TPI_GF")] <- NA

#GBDat[ !is.na( GBDat$BATHY) & GBDat$BATHY < -50 , "BATHY"] <- NA
#GBDat[ !is.na( GBDat$BATHY) & GBDat$BATHY > 0 , "BATHY"] <- NA


#GBDat_small <- GBDat[!is.na( GBDat$BATHY),]
#tmp <- colSums( GBDat_small[,c("MUZ", "SPZ", "HPZ", "NPZ")], na.rm=TRUE) 
#tmp[2] <- tmp[2] - tmp[1] # so similar amount of sites in SPZ and MUZ
#props <- tmp / nrow( GBDat_small)
#props <- props / sum( props) # 1 UP TO HERE

###################################
####  TPI to get cut points
###################################

catB <- cut(b2, breaks=Bathy.cuts, na.rm=TRUE)

plot( zones$MUZ); plot( catB, add=TRUE); plot( zones$MUZ, add=TRUE)

writeRaster(catB, "~/MBHdesignSW/design1/Bathy_cuts_SW1.tif", overwrite=TRUE)

plot(catB)


# conver to matrix for ease plotting
bathym <- raster::as.matrix(b2)
bathym
str(bathym) # 
dim(bathym) # 723 1075
bathym[70,75] # -44.75229

# transpose the axis of the matrix so x first and  y later
bathym2 <- t(bathym)
str(bathym2) # [1:1075, 1:723]
bathym2[75,70] # -44.75229

# make data frame
bathydf <- as.data.frame (
  cbind (coordinates (b2), as.numeric (bathym2)))
colnames(bathydf) <- c("Easting", "Northing", "depth")
head(bathydf)
bathydf <- bathydf[ order(bathydf$Northing, bathydf$Easting),] # order ascending first by northing and then by easting


## Setting up plotting for now and later ####
uniqueEast <- base::unique ( bathydf$Easting) # duplicate rows removed
uniqueNorth <- base::unique ( bathydf$Northing)
ELims <- range ( na.exclude ( bathydf)$Easting)
NLims <- range ( na.exclude ( bathydf)$Northing)

str(uniqueEast) ## all the x coordinates
class(uniqueEast)
str(uniqueNorth) ## all the y coordinates
str(bathym2) # the dimensions of the matrix neet to be transposed for the plot
class(bathym2)

#Fix up ordering issue
bathym2 <- bathym2[, rev ( 1 : ncol (bathym2))] # this is needed so the map looks the right way- because of the trasnposing of matrix

## plot it to see what we are dealing with ####
## these kind of plots are kind of slow ###
image.plot ( uniqueEast, uniqueNorth, bathym2,
             xlab= "Easting" , ylab= "Northing" , main= "WA South West Reefs" ,
             legend.lab= "Depth" , asp=1 , ylim= NLims, xlim= ELims,
             col=  ( tim.colors ()))     


######### Get TPI #############

#tpi <- terrain(b3, opt="TPI")
#plot(tpi)
#b3 <- raster::aggregate(b2, fact=4)
#plot(b3)
#writeRaster(tpi, "~/MBHdesignSW/spatial_data/TPIreefs12x12.tif")
tpi <- raster("~/MBHdesignSW/spatial_data/TPIreefs12x12.tif")
plot(tpi)



#Bathy.quant <- c(0,0.8,0.9,0.925,0.95,0.975,1)
tpi.quant <- c(0,0.1,0.99,1)
tpi.cuts <- quantile(tpi, tpi.quant)#c( -Inf,0.02,0.04,0.08,0.16,Inf)
#trying to make it so there is no hand-picking (except for the hand-picked function)
tmp <- cumsum( tpi.quant)
#Bathy.targetNums <- rep( floor( 18/8), 4)#floor( ( tmp / sum( tmp))[-1] * 200)#rep( 40, 5)#c( 20,20,30,65,65)
#Bathy.targetProps <-  Bathy.targetNums / sum( Bathy.targetNums)

catT <- cut( tpi, breaks=tpi.cuts, na.rm=TRUE)
plot(catT)
#writeRaster(catT,"~/MBHdesignSW/design1/TPI_cuts_SW1.tif", overwrite=TRUE)

# conver to matrix for ease plotting
tpim <- raster::as.matrix(tpi)
tpim
str(tpim) # 
dim(tpim) # 723 1075
tpim[70,75] # 0.5803415

# transpose the axis of the matrix so x first and  y later
tpim2 <- t(tpim)
str(tpim2) # [1:1075, 1:723]
tpim2[75,70] # 0.5803415

# make data frame
tpidf <- as.data.frame (
  cbind (coordinates (tpi), as.numeric (tpim2)))
colnames(tpidf ) <- c("Easting", "Northing", "tpi")
head(tpidf )
tpidf  <- tpidf [ order(tpidf $Northing, tpidf $Easting),] # order ascending first by northing and then by easting


## Setting up plotting for now and later ####
uniqueEast <- base::unique ( tpidf $Easting) # duplicate rows removed
uniqueNorth <- base::unique ( tpidf $Northing)
ELims <- range ( na.exclude ( tpidf)$Easting)
NLims <- range ( na.exclude ( tpidf )$Northing)

str(uniqueEast) ## all the x coordinates
class(uniqueEast)
str(uniqueNorth) ## all the y coordinates
str(tpim2) # the dimensions of the matrix neet to be transposed for the plot
class(tpim2)

#Fix up ordering issue
tpim2 <- tpim2[, rev ( 1 : ncol (tpim2))] # this is needed so the map looks the right way- because of the trasnposing of matrix

## plot it to see what we are dealing with ####
## these kind of plots are kind of slow ###
image.plot ( uniqueEast, uniqueNorth, tpim2,
             xlab= "Easting" , ylab= "Northing" , main= "South West Reefs" ,
             legend.lab= "TPI" , asp=1 , ylim= NLims, xlim= ELims,
             col=  ( tim.colors ()))  


#### INCLUSION PROBS ####
par ( mfrow= c ( 1 , 3 ), mar= rep ( 4 , 4 )) 

tpi.quant <- c(0,0.1,0.98,0.99,1)
tpi.cuts <- quantile(tpi, tpi.quant)#c( -Inf,0.02,0.04,0.08,0.16,Inf)
#trying to make it so there is no hand-picking (except for the hand-picked function)
tmp <- cumsum( tpi.quant)
#Bathy.targetNums <- rep( floor( 18/8), 4)#floor( ( tmp / sum( tmp))[-1] * 200)#rep( 40, 5)#c( 20,20,30,65,65)
#Bathy.targetProps <-  Bathy.targetNums / sum( Bathy.targetNums)

catT <- cut( tpi, breaks=tpi.cuts, na.rm=TRUE)
plot(catT)
#writeRaster(catT,"~/MBHdesignSW/design1/TPI_cuts_SW1.tif", overwrite=TRUE)

# set n 
n <- 36

# The number of 'depth bins' to spread sampling effort over.
nbins <- 4

# force the breaks so R doesn't use 'pretty'

breaks <- seq ( from= min ( tpidf$tpi, na.rm= TRUE ),
                to= max ( tpidf$tpi, na.rm= TRUE ), length= nbins +1 ) # -2.1618613 -0.7142973  0.7332666  2.1808306

# chech the values above with raster
minValue(tpi$tpi) # -2.161861
maxValue(tpi$tpi) # 2.180831

# Find sensible tpi bins using pre-packaged code
tmpHist <- hist ( tpidf$tpi, breaks= breaks, plot= FALSE )

#  check how it looks to see les frequent values 
# to see the plot
tmpHist <- hist ( tpidf$tpi, breaks= breaks, plot= T, freq = F )


# change breaks if needed - least interest areas should have more counts

# When reset the breaks run the tmpHist again
tmpHist <- hist ( tpidf$tpi, breaks= tpi.cuts, plot= T, freq = F )
#breaks <- c(-46.90025, -30, -20, -8.540641) #
#breaks <- c(-1.1969223, -0.8,  0.5258718,  1.9372688)
# check breaks
#tmpHist <- hist ( bathydf$depth, breaks= breaks, freq= F, plot= T )

# Find the inclusion probability for each 'stratum' (for earch 'bin')
tmpHist$inclProbs <- (n/(nbins)) / tmpHist$counts # 1.136364e-01 4.935784e-06 1.138408e-04 3.921569e-01 2.000000e+01

# Matching up locations to probabilties - in data frame
str(tpidf)
tmpHist$ID <- findInterval ( tpidf$tpi, tmpHist$breaks) # breaks coded as 1, 2, 3 and so on depending on how many bins
tmpHist$ID[20000] # 2
# not sure why the NAs, but with depth it worked fine
length(tmpHist$ID) # 48689

# see hist - quicker way to observe freq of each bin
hist(tmpHist$ID)

head(tpidf)
# A container for the design
# create data frame with each location and its inlcusion probability ####
design <- data.frame ( siteID= 1 : nrow ( tpidf),
                       Easting= tpidf$Easting, Northing= tpidf$Northing,
                       tpi= tpidf$tpi, inclProb= tmpHist$inclProbs[tmpHist$ID])
str(design)
head(design)


### test ####
str(design)
head(design)

incprob <- design[,c(2,3,5)]
head(incprob)

# make df a raster --
coordinates(incprob) <- ~Easting+Northing
gridded(incprob) <- TRUE
rasterIP <- raster(incprob)
plot(rasterIP)
cellStats(rasterIP, sum)

#writeRaster(rasterIP, "~/MBHdesignSW/design4/InclProbsDesign4n48.tif")


str(design)
length(design$Easting)
# make matrix with inclProbs from design
m <- matrix ( design$inclProb, nrow= length ( uniqueEast), byrow= F)
str(m)
head(m)
m[1]

# then plot
with ( design,
       image.plot ( uniqueEast, uniqueNorth, m,
                    xlab= "" , ylab= "" , main= "Inclusion Probability (3 clusters)" , asp= 1 ,
                    ylim= NLims, xlim= ELims))

# Take the Sample using the inclusion probabilities ####
#design$inclProb[4174928] <- 4.935784e-06 ### give this NA an inclusion value

str(design)

##### replace Nas of Inclusion probabilities for zeroes ####
names(design)
head(design)
design$inclProb[is.na(design$inclProb)] <- 0
head(design)
class(design)
any(is.na(design$inclProb))


# turn design df into matrix
designMat <- cbind(design$Easting,design$Northing, design$inclProb)
head(designMat)
str(designMat)
class(designMat)
colnames(designMat) <- c("Easting", "Northing", "inclProbs")
#Fix up ordering issue: revert ordering of columns: first column last so: inclProbs, Northing, Easting
designMat <- designMat[, rev ( 1 : ncol (designMat))]
###############

## Make data frame out of designMat to make raster ###
designDF <- as.data.frame(designMat)
head(designDF)
designDF <- designDF[ order ( designDF$Northing, designDF$Easting),] # order ascending first by northing and then by easting


# turn data frame into raster ####
# df into spatial points 
coordinates(designDF) <- ~ Easting + Northing
# coerce to SpatialPixelsDataFrame
gridded(designDF) <- TRUE
# coerce to raster
designr <- raster(designDF)
designr
plot(designr)

#second Mat - for plotting
designMat2 <- raster::as.matrix(designr, mode ='any')
dim(designMat2)
str(designMat2)
# transpose matrix
designMat3 <- t(designMat2)
dim(designMat3)
str(designMat3)
designMat3 <- designMat3[, rev ( 1 : ncol (designMat3))]

### make a new data frame out of this raster and Matrix
designdf <- as.data.frame (
  cbind ( coordinates ( designr), as.numeric ( designMat3))) 
colnames ( designdf) <- c ( "Easting" , "Northing" , "inclProbs" ) 
head(designdf)
designdf <- designdf[ order ( designdf$Northing, designdf$Easting),] # order ascending first by northing and then by easting


########## Get cluster centres #######
# Sample with 'quasiSamp' from MBH package #### this takes some time
Clusters <- quasiSamp ( n = 3, dimension= 2 ,
                    potential.sites = coordinates(designr),
                    inclusion.probs= designdf$inclProbs , nSampsToConsider= 10000) # inclProb that are not NA!
Clusters
write.csv(Clusters, "~/MBHdesignSW/design4/Clusters4.8.csv")

#### plot Clusters ####

clustersp <- Clusters
coordinates(clustersp) <- ~x+y
clustersp
proj4string(clustersp) <- proj4string(b2)

#plot(b2)
#plot(tpi)
plot(rasterIP)
plot(clustersp, pch=20, col='black',add=T)

pdf( "~/MBHdesignSW/design4/Clusters4.8_design4SW.pdf", height=7, width=8)
plot(tpi, main = "Cluster centres - SW")
points(clustersp, add=T, pch=20, col = 'black', add = TRUE)
dev.off()

###############################
####  Choose new points within clusters
####  Here I need to choose transects not points
##############################

getlocal <- function(ii){
  point <- Clusters[ii,c("x","y")]
  r2 <- rasterize( point, rasterIP, field=1)
  pbuf <- buffer( r2, width=500) ## units are in metres
  buf <- mask( rasterIP, pbuf)
  buffer <- trim(buf, pad=0)
  return( buffer)
}


sampWithOver <- 12

fullSample <- list()
fullZones <- list()

## I think in this funtion I need to change quasiSamp for TransectSamp
for( ii in 1:nrow( clustersp)){
  tmp <- getlocal(ii)
  fullZones[[ii]] <- rownames( clustersp@data)[ii]
  tmpm <- raster::as.matrix(tmp)
  tmpm <- t(tmpm)
  tmpdf <- as.data.frame (
    cbind (coordinates (tmp), as.numeric (tmpm)))
  colnames(tmpdf) <- c("x", "y", "inclProbs_design1")
  tmpdf <- tmpdf[ order(tmpdf$y, tmpdf$x),]  # order ascending first by northing and then by easting
  fullSample[[ii]] <- quasiSamp( n=sampWithOver, potential.sites=coordinates( tmp), 
                                 inclusion.probs=values( tmp), nSampsToConsider=5000)
  plot( tmp)
  points( fullSample[[ii]]$points[,c("x","y")], pch=20, col='red')
  #plot( legacySites, add=TRUE, pch=4, col='blue')
}
fullSample <- do.call( "rbind", fullSample)
fullSample$cluster <- rep( do.call( "c", fullZones), each=sampWithOver)

#fullSample$ID <- paste( fullSample$cluster, rep( paste0( "shot.",1:6), each=nrow( clustersp)), sep="_")
#fullSample2 <- SpatialPointsDataFrame( coords=fullSample[,c("x","y")], data=fullSample, proj4string=CRS(proj4string(inclProbs)))

fullSample3 <- fullSample
coordinates(fullSample3) <- ~x+y
proj4string(fullSample3) <- proj4string(tpi)

plot(rasterIP)
points(fullSample3, pch=20, col='black')
points(clustersp, pch=20, col='blue')

pdf( "~/MBHdesignSW/design3/Clusters4.8_design4SW.pdf", height=7, width=8)
plot(rasterIP)
plot(tpi)
points(fullSample3, pch=20, col='black')
points(clustersp, pch=20, col='blue')
dev.off()

### plot nicely ##
library(dichromat)
pal <- colorRampPalette(c("red","blue"))
pdf( "~/MBHdesignSW/design4/ClusteredBruvs4.8_design4SW.pdf", height=7, width=8)
plot(tpi, main ="Clustered Stereo-BRUVs - TPI", col = rev(brewer.pal(20, "RdYlBu")))
#plot(tpi, main ="Clustered Stereo-BRUVs - SW", col = pal1)
points(fullSample3, pch=20, col='black')
dev.off()

####  Write the shape files

writeOGR( fullSample3, dsn="~/MBHdesignSW/design4", layer=paste( "Bruvs4.8", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR( clustersp, dsn="~/MBHdesignSW/design4", layer=paste( "clusters4.8", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
