## Alexandra King - V00827380
## Geog418 Assignment 4

################################
#### Prepare Pollution Data ####
################################
install.packages("rgdal")
install.packages("gstat")
install.packages("sp")
install.packages("spatstat")
install.packages("maptools")
install.packages("raster")
install.packages("tmap")

library(rgdal)
library(gstat)
library(sp)
library(spatstat)  # Used for the dirichlet tessellation function
library(maptools)  # Used for conversion from SPDF to ppp
library(raster)    # Used to clip out thiessen polygons
library(tmap)

dir <- "Z:\\Geog418\\Assignment4\\Working"
setwd(dir)

#DATASET 1
#Read the pollution csv dataset.
ozone = read.csv("OZONE_PICKDATA_2016-4-30.csv", header = T, sep = ",")

#DATASET 2
#Read the monitoring station spatial dataset as an OGR data object.
monitor = readOGR(dsn = ".", layer = "airmonitoringstations")
#Extract the monitoring stations for the South Coast (SC) - taking a subset into new variable
SC.monitor = monitor[monitor$AIRBASIN %in% c("South Coast"),]
#Reproject the data to a suitable projection. Here we use a UTM projection because of the scale of the analysis.
#California is UTM Zone 11 in NAD83. (espg 26911)
SC.monitor.t = spTransform(SC.monitor, CRS("+init=epsg:26911"))


#DATASET 3
#Read the California Air Basin spatial dataset.
Ca.AirBasin = readOGR(dsn = ".", layer = "CaAirBasin")

#Extract the South Coast air basin from the spatial dataset. 
SC.AirBasin = Ca.AirBasin[Ca.AirBasin$NAME %in% c("South Coast"),] 

#Reproject the South Coast air basin spatial dataset to match the projeciton of the monitoring station dataset.  
SC.AirBasin.t = spTransform(SC.AirBasin, CRS("+init=epsg:26911"))


################################
#### Process Pollution Data ####
################################

#You need to represent each location with a single value in order to perform statistical analyses.

#Examine the first several rows of the ozone dataset. 
head(ozone)

#Looking at the date and hour columns, you can see that we need to process the data
#to get summary statistics.

#Calculate the mean and max ozone level for each site for all readings.
#only highlight and enter "aggregate(value ~ site, ozone, mean)" to get value on console
#these become just dataframes. we will need to spatialize these later on in the code
mean.ozone <- aggregate(value ~ site, ozone, mean)
max.ozone <- aggregate(value ~ site, ozone, max)

#Join the mean and max ozone values to their respective monitoring stations. In doing so, you will need to rename the 
#first column of the monitoring data to site in order to have a unique name to match the two datasets.
names(SC.monitor.t)[1] ="site"  

#Merge the the monitoring station (x) shapefile with the ozone data (y) using the site column. 
#the sp:: is saying to pull the merge tool from the sp package, instead of its default package "raster"
#'all.x = FALSE' is saying that we'll only keep the x's that match the y's, with no leftovers. but we'll keep all y's.
mrg.tab.mean <- sp::merge(SC.monitor.t, mean.ozone, by = "site", all.x = FALSE) 
mrg.tab.max <- sp::merge(SC.monitor.t, max.ozone, by = "site", all.x = FALSE)

#Create a max and a mean spatialPointDataFrame, removing N/A values
ozone.mean.spdf <- na.omit(mrg.tab.mean)
ozone.max.spdf <- na.omit(mrg.tab.max)

view(ozone.mean.spdf)

# Load and observe ozone data.... change 
tm_shape(SC.AirBasin.t) + 
  tm_polygons() +
  tm_shape(ozone.mean.spdf) +
  tm_layout(main.title = "Mean Ozone Levels in Southern Coast Air Basin") +
  tm_dots(col="value", palette = "YlOrBr", 
          title="Sampled Ozone \n(in ppm)", size=0.7) + 
  tm_legend(legend.outside=TRUE) 


#study area map
studymap_tm <- tm_shape(Ca.AirBasin) + 
    tm_fill("lightgrey") +
    tm_borders("black") +
    tm_shape(SC.AirBasin) +
    tm_fill("coral") +
    tm_borders("black") +
    tm_add_legend(type= "symbol", labels = "Southern Coast Air Basin", col="coral", shape = 19) +
    tm_layout(title = "Map of California\nAir Basins", title.position = c(0.55, 0.87),
              legend.position = c(0.55, 0.7)) +
    tm_compass(position = c(0.02, 0.085)) +
    tm_scale_bar(position= c("left", "bottom"))
studymap_tm

#SoCAB with monitoring sites
studymap2 <- tm_shape(SC.AirBasin) +
  tm_fill("coral")+
  tm_borders("black")+
  tm_add_legend


tmaptools::palette_explorer()

####################################################
### Spatial Interpolation with Thiessen Polygons ###
####################################################

# Create a tessellated surface
th  <-  as(dirichlet(as.ppp(ozone.mean.spdf)), "SpatialPolygons")

# The dirichlet function does not carry over projection information
# requiring that this information be added manually
proj4string(th) <- proj4string(ozone.mean.spdf)

# The tessellated surface does not store attribute information
# from the point data layer. We'll use the over() function (from the sp
# package) to join the point attributes to the tesselated surface via
# a spatial join. The over() function creates a dataframe that will need to
# be added to the `th` object thus creating a SpatialPolygonsDataFrame object
th.z     <- over(th, ozone.mean.spdf, fn=mean) #some get N/A values
th.spdf  <-  SpatialPolygonsDataFrame(th, th.z)

# Finally, we'll clip the tessellated  surface to the South Coast Air Basin boundaries
th.clp   <- raster::intersect(SC.AirBasin.t,th.spdf)

# Map the data
tm_shape(th.clp) + 
  tm_polygons(col="value", palette="YlOrBr",
              title="Predicted Ozone \n(in ppm)") +
  tm_layout(main.title = "Thiessen Polygons Interpolation") +
  tm_legend(legend.outside=TRUE)






########################################
#### Spatial Interpolation with IDW ####
########################################

# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(ozone.mean.spdf, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object

proj4string(grd) <- proj4string(SC.monitor.t)
P.idw <- gstat::idw(value ~ 1, ozone.mean.spdf, newdata=grd, idp=2.5)
r       <- raster(P.idw)
r.m     <- mask(r, SC.AirBasin.t)

tm_shape(r.m) + 
  tm_raster(n=10,palette = "YlOrBr",
            title="Predicted Ozone \n(in ppm)") +
  tm_layout(main.title = "IDW Interpolation") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)


#################################################
# Leave-one-out validation routine
IDW.out <- vector(length = length(ozone.mean.spdf))
for (i in 1:length(ozone.mean.spdf)) {
  IDW.out[i] <- gstat::idw(value ~ 1, ozone.mean.spdf[-i,], ozone.mean.spdf[i,], idp=2.5)$var1.pred
}

#we've run the for loop and created surfaces where we remove one point, interpolate it, and then compare
# what the interpolated value is under that point in comparison to the observation

# Plot the differences (similar to regression) using actual value of ozone to compare to interpolated IDW
#with a line of fit to see how accurate we are
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ ozone.mean.spdf$value, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ ozone.mean.spdf$value), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
#determine root mean square error
sqrt( sum((IDW.out - ozone.mean.spdf$value)^2) / length(ozone.mean.spdf))

#our plot shows that we're terrible! adjust the idp value in for loop, see how trendline and error value change

#################################################
# Implementation of a jackknife technique to estimate a confidence interval at each unsampled point.
# Create the interpolated surface
img <- gstat::idw(value~1, ozone.mean.spdf, newdata=grd, idp=2.5)
n   <- length(ozone.mean.spdf)
Zi  <- matrix(nrow = length(img$var1.pred), ncol = n)

# Remove a point then interpolate (do this n times for each point)
st <- stack()
for (i in 1:n){
  Z1 <- gstat::idw(value~1, ozone.mean.spdf[-i,], newdata=grd, idp=2.5)
  st <- addLayer(st,raster(Z1,layer=1))
  # Calculated pseudo-value Z at j
  Zi[,i] <- n * img$var1.pred - (n-1) * Z1$var1.pred
}

# Jackknife estimator of parameter Z at location j
Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )

# Compute (Zi* - Zj)^2
c1 <- apply(Zi,2,'-',Zj)            # Compute the difference
c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference

# Compute the confidence interval
CI <- sqrt( 1/(n*(n-1)) * c1)

# Create (CI / interpolated value) raster
img.sig   <- img
img.sig$v <- CI /img$var1.pred 

# Clip the confidence raster to Southern California
r <- raster(img.sig, layer="v")
r.m <- mask(r, SC.AirBasin.t)

# Plot the map
tm_shape(r.m) + tm_raster(n=7,title="95% confidence \ninterval \n(in ppm)") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)



###################################################
### Spatial Interpolation with Ordinary Krieging ##
###################################################
# just fitting on the points, not counting for trends

#assigning value = 1 for a polynomial
f.0 <- as.formula(value ~ 1) 

#pulls out the variogram, getting the mean bin values
#we can play around with our sill, range, and nugget values AND the model type
var.smpl <- variogram(f.0, ozone.mean.spdf, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=1.45e-05, model="Sph", range=40000, nugget=0))
plot(var.smpl, dat.fit, main= "Spherical Semivariogram Model: Manual")

################
# MAGIC
#here we take out the range, nugg, and sill specification ability
# we only put the model type and it automatically does the magic best fit
var.smpl <- variogram(f.0, ozone.mean.spdf, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fitEx  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model="Exp"))
plot(var.smpl, dat.fitEx, main= "Exponential Semivariogram Model: Automatic")
################


# Define the model
f.0 <- as.formula(value ~ 1) 

# Perform the krige interpolation (note the use of the variogram model
# created in the earlier step)
##ordinary kriging: not defining a polynomial trend surface, 
# it is just krieging on raw data, not looking at any trends
dat.krg <- krige( f.0, ozone.mean.spdf, grd, dat.fit)

# Convert kriged surface to a raster object for clipping
r <- raster(dat.krg)
r.m <- mask(r, SC.AirBasin.t)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=10, palette="YlOrBr",  
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_layout(main.title = "Ordinary Krieging with Spherical Model") +
  tm_legend(legend.outside=TRUE)

#can also map the variance
r   <- raster(dat.krg, layer="var1.var")
r.m <- mask(r, SC.AirBasin.t)

tm_shape(r.m) + 
  tm_raster(n=7, palette ="YlOrBr",
            title="Variance map \n(in squared ppm)") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_layout(main.title = "Variance using Spherical") +
  tm_legend(legend.outside=TRUE)

#and map the confidence interval
r   <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
r.m <- mask(r, SC.AirBasin.t)

tm_shape(r.m) + 
  tm_raster(n=7, palette ="YlOrBr",
            title="95% CI map \n(in ppm)") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_layout(main.title = "Confidence Interval using Spherical") +
  tm_legend(legend.outside=TRUE)



########################################################
### Spatial Interpolation with Trend Surface Analysis ##
########################################################

# Polynomial Trends

##############
#Define the 1st order polynomial equation
###############

#linear trend with just x + y
f.1 <- as.formula(value ~ X + Y) 

# Add X and Y to P - predict attribute based on x and y coordinates
ozone.mean.spdf$X <- coordinates(ozone.mean.spdf)[,1]
ozone.mean.spdf$Y <- coordinates(ozone.mean.spdf)[,2]

# Run the regression model (lm is a linear regression model)
# giving it value of function = x + y while looking at ozone data
lm.1 <- lm( f.1, data=ozone.mean.spdf)

# Use the regression model output to interpolate the surface
dat.1st <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.1, newdata=grd))) 

# Clip the interpolated raster to Southern California
r   <- raster(dat.1st)
r.m <- mask(r, SC.AirBasin.t)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=10, palette="YlOrBr", 
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_layout(main.title = "1st Order Polynomial Trend Surface") +
  tm_legend(legend.outside=TRUE)

#change colour pallette since this is a diverging one (do just reds or blues)
#interesting results that show ozone is increasing as we go Southeast
#ideas: going into the valleys - where wildfires occur, air is trapped in these areas
#interesting though because away from major cities.
#ocean coastal air sweeps away pollutants and dilutes the concentrations

###############
# Define the 2nd order polynomial equation
###############
#I functions means "keep this intact and don't try to interpret it" (ie there is no X variable)
f.2 <- as.formula(value ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))

# Add X and Y to P
ozone.mean.spdf$X <- coordinates(ozone.mean.spdf)[,1]
ozone.mean.spdf$Y <- coordinates(ozone.mean.spdf)[,2]

# Run the regression model again using the 2nd polynomial equation
lm.2 <- lm( f.2, data=ozone.mean.spdf)

# Use the regression model output to interpolate the surface
dat.2nd <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.2, newdata=grd))) 

# Clip the interpolated raster to South Cali
r   <- raster(dat.2nd)
r.m <- mask(r, SC.AirBasin.t)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=10, palette="YlOrBr", 
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_layout(main.title = "2nd Order Polynomial Trend Surface") +
  tm_legend(legend.outside=TRUE)

#can keep going and do third, fourth, fifth order polynomials





##################################################
## Spatial Interpolation with Universal Kriging ##
##################################################
#counting for trends in the data

#we can choose one of these trend functions to use (1st order or 2nd order) - justify your choice
f.1 <- as.formula(value ~ X + Y) 
f.2 <- as.formula(value ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))

var.smpl <- variogram(f.1, ozone.mean.spdf, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=1.45e-05, model="Sph", range=40000, nugget=0))
plot(var.smpl, dat.fit)


# Define the trend model - copy and paste the value from the equations (f1 or f2) above
f.2 <- as.formula(value ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y)) 

# Perform the krige interpolation (note the use of the variogram model
# created in the earlier step)
dat.krg <- krige( f.1, ozone.mean.spdf, grd, dat.fit)

# Convert kriged surface to a raster object for clipping
r <- raster(dat.krg)
r.m <- mask(r, SC.AirBasin.t)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=10, palette="YlOrBr",  
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_layout(main.title = "Interpolation of Ozone using Universal Krieging") +
  tm_legend(legend.outside=TRUE)

#plot the variance
r   <- raster(dat.krg, layer="var1.var")
r.m <- mask(r, SC.AirBasin.t)

tm_shape(r.m) + 
  tm_raster(n=7, palette ="YlOrBr",
            title="Variance \n(in squared ppm)") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

#plot the confidence
r   <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
r.m <- mask(r, SC.AirBasin.t)

tm_shape(r.m) + 
  tm_raster(n=7, palette ="YlOrBr",
            title="95% Confidence \nInterval \n(in ppm)") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)


