# Alexandra King - V00827380
# Geog418 - Assignment 2
# Studying Crimes in Greater Vancouver Area


#Install Libraries - don't have to do every time! (If you see "a package doesnt exist with that name"  then re-install)
install.packages("spatstat")
install.packages("rgdal")
install.packages("maptools")
install.packages("raster")
install.packages("sp")
install.packages("plyr")
install.packages("lubridate")
install.packages("tidyverse")
install.packages("grid")
install.packages("gridExtra")
install.packages("gtable")
install.packages("maps")
install.packages("bcmaps")
install.packages("tmap")
install.packages("fBasics")


#Load Libraries - if error "function not found" pops up, your libraries aren't loaded
library("spatstat")
library("rgdal")
library("maptools")
library("raster")
library("sp")
library("plyr")
library("lubridate")
library("tidyverse")
library("grid")
library("gridExtra")
library("gtable")
library("maps")
library("bcmaps")
library("tmap")
library("fBasics")


#Set Working Directory
dir <- "Z:\\Geog418\\Assignment2"
setwd(dir)
getwd()


#############################################
######## LOADING AND CLEANING DATA ########
#############################################

#Read in city shapefile
VicCity <- readOGR(dsn = ".\\Working", layer = "City_Boundary")

#here we create a variable called 'df' which is the attribute table stored in a dataframe
dfCity <- VicCity@data #Extract the attribute table as a dataframe, store in variable
class(dfCity) #ensure new variable is a dataframe
attach(dfCity) #attach dataset

View(dfCity) #allows to view attribute table


#Read in CSV
VicCrime <- read.csv(".\\Working\\Victoria_BC_Police_Department.csv")
View(VicCrime)

VicCrime$Date <- as.POSIXct(as.character(VicCrime$incident_datetime), format = "%m/%d/%Y %H:%M")
VicCrime$Year <- year(VicCrime$Date)



#clean up the columns
VicCrime_Clean <- VicCrime[,c("incident_id", "latitude", "longitude", 
                              "hour_of_day", "day_of_week", "parent_incident_type", "Date", "Year")]

VicCrime_Clean <- VicCrime_Clean[which(VicCrime_Clean$Year == 2018),]
View(VicCrime_Clean)

Coords <- VicCrime_Clean[,c("longitude", "latitude")]
crs <- CRS("+init=epsg:4326") 

VicCrimePoints <- SpatialPointsDataFrame(coords = Coords, data = VicCrime_Clean, proj4string = crs)

VicCrimePoints <- spTransform(VicCrimePoints, CRS("+init=epsg:3005"))
VicCity <- spTransform(VicCity, CRS("+init=epsg:3005"))

VicCrimePoints <- raster::intersect(VicCrimePoints, VicCity)
VicCrimePoints <- VicCrimePoints[,-c(9:14)]
View(VicCrimePoints@data)

levels(VicCrimePoints$parent_incident_type)

categories <- unique(VicCrime_Clean$parent_incident_type)
View(categories)


#create a for loop with the listName object to run through all the crimes
listNames <- as.list(unique(as.character(VicCrimePoints$parent_incident_type)))

for(i in 1:length(listNames)){
  i = 2
  kma <- VicCrimePoints[which(VicCrimePoints$parent_incident_type == listNames[[i]]),]
    
  kma$x <- coordinates(kma)[,1]
  kma$y <- coordinates(kma)[,2]
  
  #check for and remove duplicated points
  #check for duplicated points
  #finds zero distance among points
  zd <- zerodist(kma)
  zd
  
  #remove duplicates
  kma <- remove.duplicates(kma)
    
  #create an "extent" object which can be used to create the observation window for spatstat
  kma.ext <- as.matrix(extent(kma)) 
    
  #observation window
  window <- as.owin(list(xrange = kma.ext[1,], yrange = kma.ext[2,]))
    
  #create ppp oject from spatstat
  kma.ppp <- ppp(x = kma$x, y = kma$y, window = window)
    
  
  
  
  ########################################
  ######### QUADRAT ANALYSIS ##########
  ########################################
  
  ##First, determine the number of quadrats 
  quads <- 10
    
  qcount <- quadratcount(kma.ppp, nx = quads, ny = quads)
  
  png(paste("./Quad/Quadrat ", listNames[[i]], ".png",sep = ""))
  plot(kma.ppp, pch = "+", cex = 0.5, main = paste(listNames[[i]], " N = ", nrow(kma@data), sep = ""))
  plot(qcount, add = T, col = "red")
  dev.off()
  
  qcount.df <- as.data.frame(qcount)
  
  ##Second, count the number of quadrats with a distinct number of points.
  qcount.df <- plyr::count(qcount.df,'Freq')
  
  ##Change the column names so that x=number of points and f=frequency of quadrats with x point.
  colnames(qcount.df) <- c("x","f")
  
  
  #calculations for determining Variance, Mean, and VMR
  sum.f.x2 <- sum(qcount.df$f * (qcount.df$x)^2)
    
    M <- quads * quads
    
    N <- sum(qcount.df$x * qcount.df$f)
    
    sum.fx.2 <- sum((qcount.df$f * qcount.df$x)^2)
    
    
    VAR <- (sum.f.x2 - (sum.fx.2 / M)) / (M - 1)
    
    MEAN <- N / M
    
    VMR <- VAR / MEAN
    
    
    #Finally, perform the test statistic to test for the existence of a random spatial pattern.
    chi.square = VMR * (M - 1)
    p = 1 - pchisq(chi.square, (M - 1))
    
    #make a temporary data frame to store the variables in (temp because always diff results)
    #every time we loop through, a new row is added and temp values assigned to final result
    temp <- data.frame(name = listNames[[i]], VAR = VAR,
                       MEAN = MEAN, VMR = VMR, CHI = chi.square, M = M, N = N, PVAL = p)
    
    if(i == 1){
      QuadratResult <- temp
    }else{
      QuadratResult <- rbind(QuadratResult, temp)
    }
  
  
  #####################################
  ######### NEAREST NEIGHBOUR #########
  #####################################
    
  nearestNeighbour <- nndist(kma.ppp)
    
  ##Convert the nearestNeighbor object into a dataframe.
  nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
    
  ##Change the column name to "Distance"
  colnames(nearestNeighbour) = "Distance"
    
    
  ##Calculate the nearest neighbor statistic to test for a random spatial distribution.
  #mean nearest neighbour
  nnd = sum(nearestNeighbour$Distance) / nrow(nearestNeighbour)
      
    #mean nearest neighbour for random spatial distribution
      
    studyArea <- area(VicCity)
    pointDensity <- studyArea / N
      
    r.nnd = 1 / (2 * (pointDensity^-2))
      
    d.nnd = 1.07453 / (pointDensity^-2)
      
    R = nnd / r.nnd
      
    SE.NND <- 0.26136 / ((N * pointDensity)^-2)
      
    z = (nnd - r.nnd) / SE.NND
      
   
  #make another data frame table of the crime type list and store result in the NNDResult.
  temp <- data.frame(name = listNames[[i]], NND = nnd, NNDr = r.nnd, NNDd = d.nnd, R = R, Z = z)
      
      if(i == 1){
        NNDResult <- temp
      }else{
        NNDResult <- rbind(NNDResult, temp)
      }
      
  #visualize the Victoria map with points on it  
  plot(VicCity)
  plot(kma.ppp, add=TRUE)
  dev.off()
  
    
    
  #########################################
  ##########    K-FUNCTION    #############
  #########################################
  
  #basic k-function - use a K estimate function on our previously created kma.ppp object
  k.fun <- Kest(kma.ppp, correction = "Ripley")
  #plot(k.fun)
  
  #use simulation to test the point pattern against CSR - adding confidence interval bands
  k.fun.e <- envelope(kma.ppp, Kest, nsim = 99, correction = "Ripley")
  #plot(k.fun.e)
  
  png(paste("./KFunc/K_", listNames[[i]], ".png", sep = ""))
  plot(k.fun.e, main = paste(listNames[[i]], " N = ", nrow(kma@data), sep = ""))
  dev.off()
  
  
  ##########################################
  ######   KERNEL DENSITY ESTIMATION   #####
  ##########################################
  
  #2D (gaussian) kernel, compare how bandwidth (sigma) selection influences the point density estimates
  #since data are projected, sigma is represented in metres
  #eps is the width and height of the pixels (1000m X 1000m)
  #coerce to a SpatialGridDataFrame for plotting
  
  #first line does a Kernel density and puts it in a variable
  #second line puts it into a cell size spatial grid data frame
  #fourth line - takes kde.SG dataframe and adds columns for the new grid data frame from the new kde.500
  #it binds them 
  #can run cbind line multiple times with different cell sizes. **get to know cbind and rbind for future use
  
  kde.20 <- density(kma.ppp, sigma = 20, at = "pixels", eps = c(100,100))
  kde.SG <- as(kde.20, "SpatialGridDataFrame")
  kde.50 <- density(kma.ppp, sigma = 50, at = "pixels", eps = c(100,100))
  kde.SG <- cbind(kde.SG, as(kde.50, "SpatialGridDataFrame"))
  kde.100 <- density(kma.ppp, sigma = 100, at = "pixels", eps = c(100, 100))
  kde.SG <- cbind(kde.SG, as(kde.100, "SpatialGridDataFrame"))
  kde.500 <- density(kma.ppp, sigma = 500, at = "pixels", eps = c(100, 100))
  kde.SG <- cbind(kde.SG, as(kde.500, "SpatialGridDataFrame"))
  
  #sigma is the searching radius
  #cell size (1000, 1000) is the way it represents it on the screen (because KDE is a continuous surface not a grid so need a cell size to make it a grid)
  #recommend 100 by 100 for cell size of victoria. 1x1 almost killed computer
  
    
  names(kde.SG) <- c("Sig20", "Sig50", "Sig100", "Sig500")
  #plot
  #x11() #opens a new plot window
  
  png(paste("./KDE/KDE_all", listNames[[i]], ".png", sep = ""))
  spplot(kde.SG)
  dev.off()
  
  #can see how the bandwidth selection influences the density estimates
  ###can print the summary in a table like the others - we have to figure it out
  summary(kde.SG)
  
  #use cross-validation to get the bandwidth that minimizes MSE (Mean Square Error)
  bw.d <- bw.diggle(kma.ppp)
  
  #plot the "optimal" bandwidth
  png(paste("./KDE/KDE_crossVal_", listNames[[i]], ".png", sep = ""))
  plot(bw.d, ylim=c(-10, 10), main= "title")
  dev.off()
  
  #density using the cross-validation bandwidth
  kde.bwo <- as(density(kma.ppp, sigma = bw.d, at = "pixels", eps = c(100, 100)),"SpatialGridDataFrame")
  png(paste("./KDE/KDE_Best_", listNames[[i]], ".png", sep = ""))
  plot(kde.bwo, main = paste(listNames[[i]], " N = ", nrow(kma@data), sep = ""))
  dev.off()
  
}

write.csv(QuadratResult, "./Quad/QuadratResults_AllCrime_Victoria_2018.csv", row.names = FALSE)
write.csv(NNDResult, "./NND/NNDResults_AllCrime_Victoria_2018.csv", row.names = FALSE)



color=heat.colors(100, alpha=1, rev=TRUE)
hood_colors <- grey.colors(n=1, alpha=0.2)

par(oma = c(1,1,3,1))
plot(kde.bwo, col=color)
mtext("Density (points/m^2)", side=4, line=6)
plot(VicCity, col=hood_colors, add=TRUE)
title(outer=TRUE, adj=0.5, main=paste("Kernel Density Estimate of "), listNames)

