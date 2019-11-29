# Alexandra King V00827380
# Geog418 Final Project


#################################
### INSTALL AND LOAD PACKAGES ###
#################################

#install the necessary packages
install.packages("sf")
install.packages("plyr")
install.packages("dplyr")
install.packages("spdep")
install.packages("GISTools")
install.packages("raster")
install.packages("maptools")
install.packages("rgdal")
install.packages("spatstat")
install.packages("sp")
install.packages("tmap")
install.packages("gstat")
install.packages("gtable")
install.packages("grid")
install.packages("gridExtra")
install.packages("spgwr")

#Libraries
library(sf)
library(plyr)
library(dplyr)
library(spdep)
library(GISTools)
library(raster)
library(maptools)
library(rgdal)
library(spatstat)
library(sp)
library(spatstat)
library(tmap)
library(gstat)
library(gtable)
library(grid)
library(gridExtra)
library(spgwr)



###############################
### UPLOAD AND CLEAN DATA #####
###############################

# create and set working directory
dir <- "Z:\\Geog418\\FinalProject\\Working2"
setwd(dir)

# load Data 1: particulate matter dataset
pm25 <- read.csv("PM25.csv") 
# change the names of columns 1 and 2, emitting N/A values
pm25 <- pm25[,1:2]
colnames(pm25) <- c("POSTALCODE", "PM25")
pm25 <- na.omit(pm25)

# load Data 2: postal code shapefile
postalcodes <- shapefile("BC_Postal_Codes")

# load Data 3: dissemination tract and income data
income <- read.csv("Income.csv") #Read in census income data  
colnames(income) <- c("DAUID", "Income") #Select only ID and Income columns
census.tracts <- shapefile("BC_DA.shp") #Read in dissemination tract shapefile
income.tracts <- merge(census.tracts,income, by = "DAUID") #Merge income and dissemination data
nrow(income.tracts) #Determine the number of columns in the dataframe
income.tracts <- income.tracts[!is.na(income.tracts$Income),]


# Create a Classification Map of Income's distribution
# requires a shape layer, polygon information, and title
map_Income <- tm_shape(income.tracts) + 
  tm_polygons(col = "Income", 
              title = "Median Income", 
              style = "jenks", 
              palette = "YlGn", n = 8) +
  tm_layout(title = "2015 Median Incomes by Census Tract in the GVRD", title.position = c("right", "top"), 
            legend.position = c(0.015, 0.015))
map_Income


######################################
## JOINING POLLUTION DATA TO INCOME ##
######################################

# Summary:
# this section overlays the pollution points onto the census tracts and takes one 
# value per tract. This is done by first intersecting the postal codes and income tracts,
# then merging pollution data with newly updated postal code layer using codes as joining feature. 
# PM25 values are then aggregated into one value per census using the max, and the resulting
# layer is re-joined to the spatial pollution layer. The PM25 value is the max per census area.

#Select postal codes that fall within census tracts
postalcodes <- intersect(postalcodes,income.tracts)
plot(postalcodes) #See what the data looks like spatially
head(postalcodes) #See what the data looks like in tabular form
View(postalcodes@data)

#Join PM2.5 data with postal code data
pm25.spatial <- merge(postalcodes,pm25,by = "POSTALCODE")
View(pm25.spatial@data)

#Aggregate the PM2.5 values in each census in order to have a single value per DA, based on the max.
pm25.aggregate <- aggregate((as.numeric(pm25.spatial$PM25)/10)~pm25.spatial$DAUID,FUN=max)

#Re-join aggregated data to the income.tracts layer.
colnames(pm25.aggregate) <- c("DAUID", "PM25AGG") #Select only ID and Income columns
income.pm25 <- merge(income.tracts,pm25.aggregate, by = "DAUID") #Merge income and dissemination data
View(income.pm25@data)

#Re-join aggregated data to the pm25.spatial points layer.
pm25.points.aggregate <- merge(pm25.spatial, pm25.aggregate, by = "DAUID")
View(pm25.points.aggregate@data)


##################################
# Taking a Sample, Making a Grid #
##################################

#Create a subsample of the datapoints provided in the PM2.5 dataset using my provided number
set.seed(190)
sampleSize=190 
spSample <- pm25.points.aggregate[sample(1:length(pm25.points.aggregate),sampleSize),]
View(spSample@data)

# determine max and min measured values for later use in interpolation
minPM25 <- min(spSample$PM25AGG, na.rm=T)
maxPM25 <- max(spSample$PM25AGG, na.rm=T)

# Create an empty grid called grd to use in the interpolation, n is total number of cells
grd <- as.data.frame(spsample(spSample, "regular", n=40000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(spSample)



####################################################
## Step 1: Spatial Autocorrelation on Income Data ##
####################################################

# creating a neighbourhood list and a neighbourhood net using Queen's case
census.nb <- poly2nb(income.tracts)
census.net <- nb2lines(census.nb,coords=coordinates(income.tracts))

# creating a lagged weights matrix from neighbourhood list 
census.lw <- nb2listw(census.nb, zero.policy = TRUE, style = "W")
print.listw(census.lw, zero.policy = TRUE)


################################
# GLOBAL MORAN'S I 

# uses variable Income and neighbourhood weights matrix list
mi <- moran.test(income.tracts$Income, census.lw, zero.policy = TRUE)
mi

# determining the range of morans I
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(census.lw)

# obtaining values for mI, eI, and var from the matrix
mI <- mi$estimate[[1]]
eI <- mi$estimate[[2]]
var <- mi$estimate[[3]]
# calculating the z-score
z <- (mI - eI) / sqrt(var)

# determine significance using 95% confidence
if(z > 1.96 | z < -1.96){
  sig <- "Yes"
}else{
  sig <- "No"
}

# determine spatial pattern based on mI value and significance
if(sig == "Yes"){
  if(mI > eI){
    pattern <- "Clustered"
  }else{
    pattern <- "Dispersed"
  }
}else{
  pattern <- "Random"
}


##############################
# LOCAL MORAN's I 

# run a local Moran's I test
lisa.test <- localmoran(income.tracts$Income, census.lw)

# store results in new columns in the dataset
income.tracts$Ii <- lisa.test[,1]
income.tracts$E.Ii<- lisa.test[,2]
income.tracts$Var.Ii<- lisa.test[,3]
income.tracts$Z.Ii<- lisa.test[,4]
income.tracts$P<- lisa.test[,5]

# determine Moran's I legend parameters
minIi <- min(income.tracts$Ii, na.rm=T)
maxIi <- max(income.tracts$Ii, na.rm=T)
expIi <- lisa.test[1,2]

# create map of local Moran's I results
map_LISA <- tm_shape(income.tracts) + 
  tm_polygons(col = "Ii", 
              title = "Local Moran's I", 
              style = "fixed",
              breaks = c(minIi, expIi-0.005, expIi+0.005, maxIi),
              labels = c("Negative Autocorrelation", "Random",
                         "Positive Autocorrelation"),
              midpoint = NA,
              palette = "BrBG", n = 3) +
  tm_layout(title = "Local Moran's I for Income in Vancouver Census Tracts", 
            legend.position = c(0.015, 0.015))
map_LISA


# map the p-values to visualize significance
# correct p-values for a two-tail significance test
income.tracts$P2Tail <- income.tracts$P * 2

# determine p-value legend parameters
minP <- min(income.tracts$P2Tail, na.rm=T)
maxP <- max(income.tracts$P2Tail, na.rm=T)

#create map
map_PVals <- tm_shape(income.tracts) + 
  tm_polygons(col = "P2Tail", 
              title = "P Values", 
              style = "fixed",
              breaks = c(minP, 0.05, maxP),
              labels = c("p < 0.05", "p > 0.05"),
              palette = "RdGy", n = 3) +
  tm_layout(title = "P values from Income Local Moran's I Test", 
            legend.position = c(0.015, 0.015))
map_PVals


# graph Local Moran's I results in a scatterplot
plot.Inc <- moran.plot(income.tracts$Income, census.lw, 
                       zero.policy=NULL, spChk=NULL, labels=NULL, 
                       main= "Scatterplot of Income Spatial Autocorrelation by Census Tracts", 
                       xlab="Income", ylab="Spatially Lagged Income", quiet=NULL)




#######################################################
## Step 2: Spatial Interpolation on Pollution Points ##
#######################################################

#########################
# IDW INTERPOLATION

P.idw <- gstat::idw(PM25AGG ~ 1, spSample, newdata=grd, idp=2.5)
r       <- raster(P.idw)
r.m     <- mask(r, income.tracts)


# map the IDW results
tm_shape(r.m) + 
  tm_raster(n=10,palette = "Reds",
            title="Predicted PM2.5 \n(in ug/m3)") +
  tm_layout(main.title = "IDW Interpolation") +
  tm_shape(spSample) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)


# Leave-one-out validation routine
IDW.out <- vector(length = length(spSample))
for (i in 1:length(spSample)) {
  IDW.out[i] <- gstat::idw(PM25AGG ~ 1, spSample[-i,], spSample[i,], idp=6)$var1.pred
}


# Plot the differences (similar to regression) using actual value of pm2.5 to compare to interpolated IDW
#with a line of fit to see how accurate we are
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ spSample$PM25AGG, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ spSample$PM25AGG), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
#determine root mean square error
sqrt( sum((IDW.out - spSample$PM25AGG)^2) / length(spSample))


# Implementation of a jackknife technique to estimate a confidence interval at each unsampled point.
# Create the interpolated surface
img <- gstat::idw(PM25AGG~1, spSample, newdata=grd, idp=2.5)
n   <- length(spSample)
Zi  <- matrix(nrow = length(img$var1.pred), ncol = n)

# Remove a point then interpolate (do this n times for each point)
st <- stack()
for (i in 1:n){
  Z1 <- gstat::idw(PM25AGG~1, spSample[-i,], newdata=grd, idp=2.5)
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

# Clip the confidence raster to GVRD
r <- raster(img.sig, layer="v")
r.m <- mask(r, income.tracts)

# Plot the map
tm_shape(r.m) + tm_raster(n=7,title="95% confidence \ninterval \n(in ppm)") +
  tm_shape(spSample) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)



#################################
# TREND SURFACE ANALYSIS

#Define the 1st order polynomial equation
f.1 <- as.formula(PM25AGG ~ X + Y) 

# Add X and Y to P - predict attribute based on x and y coordinates
spSample$X <- coordinates(spSample)[,1]
spSample$Y <- coordinates(spSample)[,2]

# Run the regression model (lm is a linear regression model)
# giving it value of function = x + y while looking at PM2.5 data
lm.1 <- lm( f.1, data=spSample)

# Use the regression model output to interpolate the surface
dat.1st <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.1, newdata=grd))) 

# Clip the interpolated raster to Vancouver
r   <- raster(dat.1st)
r.m <- mask(r, income.tracts)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=10, palette="YlOrBr", 
            title="Predicted PM2.5 \n(ug/m3)") +
  tm_shape(spSample) + tm_dots(size=0.15) +
  tm_layout(main.title = "1st Order Polynomial Trend Surface") +
  tm_legend(legend.outside=TRUE)


#Define the 2nd order polynomial equation
f.2 <- as.formula(PM25AGG ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))

# Add X and Y to P
spSample$X <- coordinates(spSample)[,1]
spSample$Y <- coordinates(spSample)[,2]

# Run the regression model again using the 2nd polynomial equation
lm.2 <- lm( f.2, data=spSample)

# Use the regression model output to interpolate the surface
dat.2nd <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.2, newdata=grd))) 

# Clip the interpolated raster to GVRD income tracts
r   <- raster(dat.2nd)
r.m <- mask(r, income.tracts)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=10, palette="YlOrBr", 
            title="Predicted PM2.5 \n(ug/m3)") +
  tm_shape(spSample) + tm_dots(size=0.15) +
  tm_layout(main.title = "2nd Order Polynomial Trend Surface") +
  tm_legend(legend.outside=TRUE)




########################################
# UNIVERSAL KRIGING 

# choosing to use the 2nd order function developed above
f.2 <- as.formula(PM25AGG ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))

# choosing to set my own range, sill, and nugget for my model of best fit
var.smpl <- variogram(f.2, spSample, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=0.93, model="Exp", range=5, nugget=0))
plot(var.smpl, dat.fit)

# Define the trend model again
f.2 <- as.formula(PM25AGG ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y)) 

# Perform the krige interpolation using trend function and variogram created above
dat.krg <- krige( f.2, spSample, grd, dat.fit)
View(dat.krg@data)

# rename the interpolated values from var1.pred to PM25
colnames(dat.krg@data)[1] <- "PM25"
View(dat.krg@data)

# Convert kriged surface to a raster object for clipping
r <- raster(dat.krg)
r.m <- mask(r, income.tracts)

# Plot the map showing interpolation results
tm_shape(r.m) + 
  tm_raster(n=10, palette="YlOrBr",  
            title="Predicted PM2.5 \n(ug/m3)") +
  tm_shape(spSample) + tm_dots(size=0.15) +
  tm_layout(main.title = "Interpolation of PM2.5 using Universal Krieging") +
  tm_legend(legend.outside=TRUE)


# plot the variance
rvar   <- raster(dat.krg, layer="var1.var")
r.m.var <- mask(rvar, income.tracts)

tm_shape(r.m.var) + 
  tm_raster(n=7, palette ="YlOrBr",
            title="Variance \n(ug/m3)") +tm_shape(spSample) + tm_dots(size=0.15) +
  tm_legend(legend.outside=TRUE)

# plot the confidence
rconf   <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
r.m.conf <- mask(rconf, income.tracts)

tm_shape(r.m.conf) + 
  tm_raster(n=7, palette ="YlOrBr",
            title="95% Confidence \nInterval \n(ug/m3)") +tm_shape(spSample) + tm_dots(size=0.15) +
  tm_legend(legend.outside=TRUE)




###############################
## Step 3: LINEAR REGRESSION ##
###############################

##############################
# Data Prep - combining outputs from interpolation with income data to prep for regression

# could reduce number of cells by aggregating values (fact = ), but chose not to
step.1 <- aggregate(r.m, fact=1, fun=mean)
plot(step.1)

# convert the raster dataset to points
step.2 <-  rasterToPoints(step.1,fun=NULL, spatial=FALSE, crs=census.tracts)

# convert the point dataset to a spatial dataframe
step.2 <- as.data.frame(step.2) 

# assign coordinates to a new object, utilizing an existing projection from census.tracts
Coords <- step.2[,c("x", "y")] 
crs <- crs(census.tracts) 

# create a spatial points dataframe
step.3 <- SpatialPointsDataFrame(coords = Coords, data = step.2, proj4string = crs) 

# aggregate points into census tracts
step.4 <- aggregate(x=step.3,by=income.tracts, FUN=mean) 

#get the intersection of step.4 with the income.tracts dataset (this will take a while)
step.5 <- intersect(step.4,income.tracts)   
pm.income.poly <- step.5

#remove the values that were interpolated outside of the 'measured points' range, and N/A values
pm.income.poly <- pm.income.poly[!is.na(pm.income.poly$PM25),]
pm.income.poly <- subset(pm.income.poly, pm.income.poly$PM25 >= minPM25 & pm.income.poly$PM25 <= maxPM25)
View(pm.income.poly@data)


################################
## Linear Regression 

#Plot income and PM2.5 from the pm.income.poly dataset you created
plot(pm.income.poly$PM25~pm.income.poly$Income)

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
pm.income.poly <-  pm.income.poly[pm.income.poly$PM25 != 0, ]

# plot the data again
plot(pm.income.poly$PM25~pm.income.poly$Income)

#Perform a linear regression on the two variables, using PM2.5 as dependent
lm.model <- lm(pm.income.poly$PM25~pm.income.poly$Income)

#Add the regression model to the plot you created
abline(lm.model)

#Get the summary of the results
summary(lm.model)


#################################
# Acquiring Residuals

# obtain the residuals from the model
model.resids <- as.data.frame(residuals.lm(lm.model))

# add the residuals to your spatialpolygon dataframe
pm.income.poly$residuals <- residuals.lm(lm.model)

# observe the result to make sure it looks correct
head(pm.income.poly)




###################################################
## Step 4: Spatial Autocorrelation of Residuals ##
###################################################

# creating a neighbourhood list and a neighbourhood net using Queen's case
pmincome.nb <- poly2nb(pm.income.poly)
pmincome.net <- nb2lines(pmincome.nb,coords=coordinates(census.tracts))

# creating a lagged weights matrix from neighbourhood list 
pmincome.lw <- nb2listw(pmincome.nb, zero.policy = TRUE, style = "W")
print.listw(pmincome.lw, zero.policy = TRUE)


#########################
# GLOBAL MORAN'S I 

# uses variable Residuals and neighbourhood weights matrix list
miRes <- moran.test(pm.income.poly$residuals, pmincome.lw, zero.policy = TRUE)
miRes

# determining the range of morans I
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(pmincome.lw)

# obtaining values for mI, eI, and var from the matrix
mIres <- miRes$estimate[[1]]
eIres <- miRes$estimate[[2]]
varRes <- miRes$estimate[[3]]
# calculating the z-score
zRes <- (mIres - eIres) / sqrt(varRes)

# determine significance using 95% confidence
if(zRes > 1.96 | zRes < -1.96){
  sigRes <- "Yes"
}else{
  sigRes <- "No"
}

# determine spatial pattern based on mI value and significance
if(sigRes == "Yes"){
  if(mIres > eIres){
    patternRes <- "Clustered"
  }else{
    patternRes <- "Dispersed"
  }
}else{
  patternRes <- "Random"
}


##########################
# LOCAL MORAN's I 

# run a local Moran's I test
lisa.test.res <- localmoran(pm.income.poly$residuals, pmincome.lw)

# store results in new columns in the dataset
pm.income.poly$Ii.res <- lisa.test.res[,1]
pm.income.poly$E.Ii.res<- lisa.test.res[,2]
pm.income.poly$Var.Ii.res<- lisa.test.res[,3]
pm.income.poly$Z.Ii.res<- lisa.test.res[,4]
pm.income.poly$P.res<- lisa.test.res[,5]

View(pm.income.poly@data)

# map the local Moran's i results
# determine Moran's I legend parameters
minIi.res <- min(pm.income.poly$Ii.res, na.rm=T)
maxIi.res <- max(pm.income.poly$Ii.res, na.rm=T)
expIi.res <- lisa.test.res[1,2]

# create map
map_LISAres <- tm_shape(pm.income.poly) + 
  tm_polygons(col = "Ii", 
              title = "Local Moran's I", 
              style = "fixed",
              breaks = c(minIi.res, expIi.res-0.005, expIi.res+0.005, maxIi.res),
              labels = c("Negative Autocorrelation", "Random",
                         "Positive Autocorrelation"),
              midpoint = NA,
              palette = "BrBG", n = 3) +
  tm_layout(title = "Local Moran's I of the Residuals", 
            legend.position = c(0.015, 0.015))
map_LISAres


# map the p-values to visualize significance
# correct p-values for a two-tail significance test
pm.income.poly$P2Tail.res <- pm.income.poly$P.res * 2

# determine p-value legend parameters
minP.res <- min(pm.income.poly$P2Tail.res, na.rm=T)
maxP.res <- max(pm.income.poly$P2Tail.res, na.rm=T)

#create map
map_PValsRes <- tm_shape(pm.income.poly) + 
  tm_polygons(col = "P2Tail", 
              title = "P Values", 
              style = "fixed",
              breaks = c(minP.res, 0.05, maxP.res),
              labels = c("p < 0.05", "p > 0.05"),
              palette = "RdGy", n = 3) +
  tm_layout(title = "P values from Local Moran's I Test on Residuals", 
            legend.position = c(0.015, 0.015))
map_PValsRes


# graph Local Moran's I results in a scatterplot
plot.Residuals <- moran.plot(pm.income.poly$residuals, pmincome.lw, 
                             zero.policy=NULL, spChk=NULL, labels=NULL, 
                             main= "Scatterplot of Residuals Autocorrelation", 
                             xlab="Residuals", ylab="Spatially Lagged Residuals", quiet=NULL)




#################################################
## Step 5: Geographically Weighted Regression ###
#################################################

# add the polygon coordinates to the spatialpolygondataframe.
# obtain the polygon coordinates using the "coordinates" function from the sp library
pm.income.poly.coords <- sp::coordinates(pm.income.poly)

# observe the result
head(pm.income.poly.coords)

# add the polygon coordinates to the spatialpolygondataframe
pm.income.poly$X <- pm.income.poly.coords[,1]
pm.income.poly$Y <- pm.income.poly.coords[,2]
head(pm.income.poly)

# determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(pm.income.poly$PM25~pm.income.poly$Income, 
                        data=pm.income.poly, coords=cbind(pm.income.poly$X,pm.income.poly$Y),adapt=T) 

# perform GWR on the two variables with the bandwidth determined above (takes a looong while)
gwr.model = gwr(pm.income.poly$PM25~pm.income.poly$Income, 
                data=pm.income.poly, coords=cbind(pm.income.poly$X,pm.income.poly$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

# print the results of the model
gwr.model

# observe the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)


# add the local r2 values (removing negatives) and coefficients to the original layer
pm.income.poly$localr <- results$localR2
pm.income.poly$coeff <- results$pm.income.poly.Income
pm.income.poly <- subset(pm.income.poly, pm.income.poly$localr >= 0)

View(pm.income.poly@data)

#tmap_mode("view") 

# create map of r-squared values
map_Localr <- tm_shape(pm.income.poly) + 
  tm_polygons(col = "localr", 
              title = "Local R2 Values", 
              style = "jenks", 
              palette = "Blues", n = 8) +
  tm_layout(title = "Local R-squared Values", 
            legend.position = c(0.015, 0.015)) 
map_Localr



# create map of coefficients
map_Coeff <- tm_shape(pm.income.poly) + 
  tm_polygons(col = "coeff", 
              title = "Coefficients", 
              style = "jenks", 
              palette = "Oranges", n = 8) +
  tm_layout(title = "Coefficient", 
            legend.position = c(0.015, 0.015)) 
map_Coeff


View(pm.income.poly@data)
View(results)



####################################
## Step 6: Point Pattern Analysis ##
####################################

######################################
## QUADRAT ANALYSIS 

# re-project pollution sample dataset and income tracts to BC Albers
spSample.t <- spTransform(spSample, CRS("+init=epsg:3005"))
income.tracts.t <- spTransform(income.tracts, CRS("+init=epsg:3005"))

# add coordinates
spSample.t$x <- coordinates(spSample.t)[,1]
spSample.t$y <- coordinates(spSample.t)[,2]
View(spSample.t@data)

# check for and remove duplicated points
zd <- zerodist(spSample.t)
zd
spSample.t <- remove.duplicates(spSample.t)

# create an "extent" object which can be used to create the observation window for spatstat
spSample.ext <- as.matrix(extent(spSample.t)) 

# create an observation window
window <- as.owin(list(xrange = spSample.ext[1,], yrange = spSample.ext[2,]))

# create a ppp object from spatstat
spSample.ppp <- ppp(x = spSample.t$x, y = spSample.t$y, window = window)

# determine the number of quadrats for the analysis
quads <- 15
qcount <- quadratcount(spSample.ppp, nx = quads, ny = quads)

# plot the output
plot(spSample.ppp, pch = "+", cex = 0.5, main = "Quadrat Analysis on PM2.5 Points")
plot(qcount, add = T, col = "red")
plot(income.tracts.t, add = TRUE)

qcount.df <- as.data.frame(qcount)

##Second, count the number of quadrats with a distinct number of points.
qcount.df <- plyr::count(qcount.df,'Freq')

##Change the column names so that x=number of points and f=frequency of quadrats with x point.
colnames(qcount.df) <- c("x","f")


# calculations for determining Variance, Mean, and VMR
sum.f.x2 <- sum(qcount.df$f * (qcount.df$x)^2)

M <- quads * quads

N <- sum(qcount.df$x * qcount.df$f)

sum.fx.2 <- sum((qcount.df$f * qcount.df$x)^2)

QuadVAR <- (sum.f.x2 - (sum.fx.2 / M)) / (M - 1)

QuadMEAN <- N / M

QuadVMR <- QuadVAR / QuadMEAN


# perform a chi-square test to test for the significant existence of a random spatial pattern.
chi.square = QuadVMR * (M - 1)
p = 1 - pchisq(chi.square, (M - 1))



################################################
###KERNEL DENSITY ESTIMATION

# 2D (gaussian) kernel, compare how bandwidth (sigma) selection influences the point density estimates
# since data are projected, sigma is represented in metres
# eps is the width and height of the pixels (100m X 100m)
# coerce to a SpatialGridDataFrame for plotting

kde.100 <- density(spSample.ppp, sigma = 100, at = "pixels", eps = c(100, 100))
kde.SG <- as(kde.100, "SpatialGridDataFrame")
kde.500 <- density(spSample.ppp, sigma = 500, at = "pixels", eps = c(100, 100))
kde.SG <- cbind(kde.SG, as(kde.500, "SpatialGridDataFrame"))

# can see how the bandwidth selection influences the density estimates
summary(kde.SG)

# use cross-validation to get the bandwidth that minimizes MSE
bw.d <- bw.diggle(spSample.ppp)

# plot the "optimal" bandwidth
plot(bw.d, ylim=c(-10, 10), main= "title")

# plot the density using the cross-validation bandwidth
kde.bwo <- density(spSample.ppp, sigma = bw.d, at = "pixels", eps = c(100, 100))
plot(kde.bwo, main = "Kernel Density Estimation on PM2.5 Points")

# visualize the GVRD map with PM2.5 points on it to compare to KDE surface
plot(income.tracts.t)
plot(spSample.ppp, add=TRUE)
dev.off()

