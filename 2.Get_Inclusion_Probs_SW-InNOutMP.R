# GEt inclusion probabilities ---

## Create inclusion probabilities ####

library( rgdal)
library( sp)
library( raster)

# clear environment ----
rm(list = ls())

# Set working directory ####
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
s.dir <- paste(w.dir, "spatial_data", sep ='/')
d.dir <- paste(w.dir, "InNOutMP", sep ='/')



#read in data

SWDat <- readRDS( "SWData_forInNOutMP.RDS")
sw_rasters <- readRDS("SWRasters_forInNOutMP.RDS")
zones <- readRDS( "SWZones_forInNOutMP.RDS")

# SW marine park polygon ----
swnp <- readOGR(paste(s.dir, "SW_CMR_NP.shp", sep='/'))
plot(swnp, add=T)

####  Straw man for numbers of samples in each region ----

straw.nums <- c(40, 24)  #numbers of drops in and out
straw.props <- straw.nums / sum( straw.nums) # 0.625 0.375
names( straw.nums) <- names( straw.props) <- c("In", "Out")
saveRDS( straw.nums, file="StrawmanNumbers_Zones.RDS")



####  Hand-picking Bathy cut points ----
####  And their numbers of drops

#Bathy.quant <- c(0,0.8,0.9,0.925,0.95,0.975,1)
#Bathy.quant <- c(0,0.25,0.5,0.8,0.9,0.925,0.95,0.975,1)
Bathy.quant <- c(0,0.009,0.06,0.07,0.48,0.5,0.8,0.85,1)
Bathy.cuts <- quantile( sw_rasters$bathy, Bathy.quant)#c( -Inf,0.02,0.04,0.08,0.16,Inf)
#trying to make it so there is no hand-picking (except for the hand-picked function)
tmp <- cumsum( Bathy.quant)
#Bathy.targetNums <- rep( floor( 18/8), 4)#floor( ( tmp / sum( tmp))[-1] * 200)#rep( 40, 5)#c( 20,20,30,65,65)
Bathy.targetNums <- rep( floor( 64/8), 2) # 8 8
Bathy.targetProps <-  Bathy.targetNums / sum( Bathy.targetNums) # 0.5 0.5


# Proportion of potential sites in each zone ----

SWDat_small <- SWDat[!is.na( SWDat$Bathy),]
tmp <- colSums( SWDat_small[,c("InsideMP", "OutsideMP")], na.rm=TRUE) 
#tmp[2] <- tmp[2] - tmp[1] # so similar amount of sites in SPZ and MUZ
tmp[1] # 7428 
tmp[2] # 2623 
props <- tmp / nrow( SWDat_small) # inside 0.6966145 - outside 0.2459908 
props <- props / sum( props) # inside 0.7390309 - outside 0.2609691 


###################################
####  To get cut points
###################################

catB <- cut( sw_rasters$bathy, breaks=Bathy.cuts, na.rm=TRUE)

plot( zones$InsideMP, add=T); plot( catB, add=TRUE); plot( zones$OutsideMP, add=TRUE)
plot(swnp, add=T)

writeRaster(catB, file='Bathy_cuts_InNOutMP-d1.tif', overwrite=TRUE)

plot(catB)


##################################
####  Within each zone (incl probs)
####  Weight according to straw.props
##################################



inclProbs <- catB
for( zz in c( "InsideMP", "OutsideMP")){
  print( zz)
  #if( zz == "MUZ")
  #zoneID <- extract( x=catB, y=zones$MUZ, cellnumbers=TRUE)
  #zoneID <- extract( x=catB, y=zones$MUZ-zones$NPZ, cellnumbers=TRUE)
  #else
  zoneID <- extract( x=catB, y=zones[[zz]], cellnumbers=TRUE)
  propsOfbathy <- table( catB@data@values[zoneID[[1]][,"cell"]])
  propsOfbathy <- propsOfbathy / sum( propsOfbathy)
  tmp <- Bathy.targetProps / propsOfbathy #the desired inclusion probs (unstandardised)
  for( ii in 1:length( propsOfbathy)){
    inclProbs[zoneID[[1]][,"cell"]][zoneID[[1]][,"value"]==ii] <- tmp[ii]
  }
  inclProbs[zoneID[[1]][,"cell"]][is.na( inclProbs[zoneID[[1]][,"cell"]])] <- 0
  inclProbs[zoneID[[1]][,"cell"]] <- inclProbs[zoneID[[1]][,"cell"]] / sum( inclProbs[zoneID[[1]][,"cell"]])
}
inclProbs@data@values[inclProbs@data@values %in% c(0,1,2,3,4,5,6,7,8)] <- NA  #cheats way to crop
plot( inclProbs)


#standardising so that the zone totals are correct according to straw.props | straw.nums
InMP <- extract( x=catB, y=zones$InsideMP, cellnumbers=TRUE)
OutMP <- extract( x=catB, y=zones$OutsideMP, cellnumbers=TRUE)

#HPZZone <- extract( x=catB, y=zones$HPZ, cellnumbers=TRUE)
#NPZZone <- extract( x=catB, y=zones$NPZ, cellnumbers=TRUE)

inclProbs@data@values[InMP[[1]][,'cell']] <- inclProbs@data@values[InMP[[1]][,'cell']] * straw.props["In"]
inclProbs@data@values[OutMP[[1]][,'cell']] <- inclProbs@data@values[OutMP[[1]][,'cell']] * straw.props["Out"]

#inclProbs@data@values[HPZZone[[1]][,'cell']] <- inclProbs@data@values[HPZZone[[1]][,'cell']] * straw.props["HPZ"]
#inclProbs@data@values[NPZZone[[1]][,'cell']] <- inclProbs@data@values[NPZZone[[1]][,'cell']] * straw.props["NPZ"]

plot(inclProbs)


writeRaster( inclProbs, file='inclProbs_InNOutMP-d1.tif', overwrite=TRUE)
