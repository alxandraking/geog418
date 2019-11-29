### Alexandra King - V00827380
### GEOG 418 - LAB 3 - SPATIAL AUTOCORRELATION ###


#################################
### INSTALL AND LOAD PACKAGES ###
#################################

#install the necessary packages
install.packages("plyr")
install.packages("dplyr")
install.packages("spdep")
install.packages("GISTools")
install.packages("raster")
install.packages("maptools")
install.packages("rgdal")
install.packages("tmap")
install.packages("BAMMtools")
install.packages("shinyjs")
install.packages("gtable")
install.packages("grid")
install.packages("gridExtra")

#load the packages
library(plyr)
library(dplyr)
library(spdep)
library(GISTools)
library(raster)
library(maptools)
library(rgdal)
library(tmap)
library(BAMMtools)
library(shinyjs)
library(gtable)
library(grid)
library(gridExtra)


###############################
### UPLOAD AND CLEAN DATA #####
###############################

# create and set your working directory
dir <- "Z:\\Geog418\\Assignment3"
  setwd(dir)


# upload two datasets, one shapefile for tract boundaries, 
# and one csv file with census info
tracts <- readOGR(dsn = ".\\Working", layer = "Vic_Census")
census.16 <- read.csv(".\\Working\\CensusTractData.csv")


# merge the two datasets together based on a common field in each
# in this case, GUID is equal in each dataset
crd.data <- merge(tracts, census.16, by = "GUID")


# remove any items in your columns of interest with an N/A value
crd.data <- crd.data[!is.na(crd.data$MdInc),]
crd.data <- crd.data[!is.na(crd.data$RefIm),]

# classify and view a summary of you new layer
class(crd.data)
summary(crd.data)



#######################################
## CLASSIFICATION MAP OF A VARIABLE ###
#######################################

# a useful function to explore colour pallette options
tmaptools::palette_explorer()


# creating a map to display your variable's distribution
# requires a shape layer, polygon information, and title
map_MedInc <- tm_shape(crd.data) + 
  tm_polygons(col = "MdInc", 
              title = "Median Income", 
              style = "jenks", 
              palette = "YlGn", n = 6) +
  tm_layout(title = "Median Income Distribution in the Victoria CRD in 2015")


# generate the map
map_MedInc






########################


map_Refugee <- tm_shape(crd.data) + 
  tm_polygons(col = "RefIm", 
              title = "Refugee Population", 
              style = "jenks", 
              palette = "YlGn", n = 6) +
  tm_layout(title = "Refugee Populations in the Victoria CRD in 2015")
map_Refugee



######################################
#### NEIGHBOURHOOD WEIGHTS MATRIX ####
######################################


# creating a neighbourhood list and a neighbourhood net
# Queen's case
crd.nbQueen <- poly2nb(crd.data)
crd.netQ <- nb2lines(crd.nbQueen,coords=coordinates(crd.data))
# Rook's case
crd.nbRook <- poly2nb(crd.data, queen = FALSE)
crd.netR <- nb2lines(crd.nbRook,coords=coordinates(crd.data))


# mapping an individual neighbourhood net
tm_shape(crd.data) + tm_borders(col='darkgrey') + 
  tm_shape(crd.netQ) + tm_lines(col='red', lwd = 2) +
  tm_layout(title = "Neighbourhood Net using Queen's Case")

# mapping both together
tm_shape(crd.data) + tm_borders(col='darkgrey') + 
  tm_shape(crd.netQ) + tm_lines(col='red', lwd = 2) +
  tm_shape(crd.netR) + tm_lines(col='lightblue', lwd = 2) +
  tm_add_legend(type = "line", labels = c("Queens", "Rooks"), 
                col = c(adjustcolor("red", alpha.f = 1),adjustcolor("lightblue", alpha.f = 1)), 
                lwd = c(5,5)) +
  tm_layout(title = "Neighbourhood Net Comparing Queen and Rook")



# creating a lagged weights matrix from neighbourhood list 
crd.lw <- nb2listw(crd.nbRook, zero.policy = TRUE, style = "W")
print.listw(crd.lw, zero.policy = TRUE)







crd.lw2 <- nb2listw(crd.nbQueen, zero.policy = TRUE, style = "W")
print.listw(crd.lw2, zero.policy = TRUE)

crd.lwU <- nb2listw(crd.nbRook, zero.policy = TRUE, style = "U")
print.listw(crd.lw, zero.policy = TRUE)

##############################
#### MAPPING LAGGED MEANS ####
##############################

# calculate lagged means using the weights matrix list and attribute's column
crd.data$IncLagMeans = lag.listw(crd.lw, crd.data$MdInc, zero.policy = TRUE)

# determine difference between individuals and their neighbourhoods
crd.data$difLagMeans = abs((crd.data$IncLagMeans - crd.data$MdInc))

# mapping the difference in lagged means
map_LagMean <- tm_shape(crd.data) + 
  tm_polygons(col = "difLagMeans", 
              title = "Median Income\nLagged Means",
              style = "jenks", 
              palette = "Oranges", n = 8) +
  tm_layout(title = "Lagged Means for Median Income Distribution")
map_LagMean



#REFUGEES
#calculate the lag means using the weighted matrix and our data
crd.data$RefLagMeans = lag.listw(crd.lw, crd.data$RefIm, zero.policy = TRUE)
crd.data$difRefLagMeans = abs((crd.data$RefLagMeans - crd.data$RefIm))

map_LagMean2 <- tm_shape(crd.data) + 
  tm_polygons(col = "difRefLagMeans", 
              title = "Refugee Immigration\nLagged Means", 
              style = "jenks", 
              palette = "Oranges", n = 8) +
  tm_layout(title = "Lagged Means for Refugee Populations")


map_LagMean2


########################
### GLOBAL MORANS I ####
########################

# running the global moran's I statistic using the variable data
# and weights matrix neighbourhood list
mi <- moran.test(crd.data$MdInc, crd.lw, zero.policy = TRUE)
mi


# determining the range of morans I
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(crd.lw)


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




#####
#REFUGEES

mi2 <- moran.test(crd.data$RefIm, crd.lw, zero.policy = TRUE)
mi2

#getting the full range of morans i
#making function and passing it the lagged weight matrix (lw)
#looking at neighbours and producing a range of expected moran's i values - what they could be
#so you can then interpret what you're going to get. gives you the full range of possibilities
moran.range2 <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range2(crd.lw)

# looking at variables individually - taking these from the output matrix 
# generated above
mI2 <- mi2$estimate[[1]]
eI2 <- mi2$estimate[[2]]
var2 <- mi2$estimate[[3]]

# write your z-score formula (on slides)
z2 <- (mI2 - eI2) / sqrt(var2)



# determine significance through a conditional statement
# based on z-score with a 95% confidence, two-tailed
if(z > 1.96 | z < -1.96){
  sig <- "Yes"
}else{
  sig <- "No"
}

if(z2 > 1.96 | z < -1.96){
  sig2 <- "Yes"
}else{
  sig2 <- "No"
}

# determine cluster or dispersed pattern based on morans I
if(sig == "Yes"){
  if(mI > eI){
    pattern <- "Clustered"
  }else{
    pattern <- "Dispersed"
  }
}else{
    pattern <- "Random"
}

if(sig2 == "Yes"){
  if(mI2 > eI2){
    pattern2 <- "Clustered"
  }else{
    pattern2 <- "Dispersed"
  }
}else{
  pattern2 <- "Random"
}



###########################
##### CREATE A TABLE ######
###########################

# create objects for the table columns
variables = c("Median Income", "Refugee Population")
moransI = c(mI, mI2) 
expectedI = c(eI, eI2) 
variance = c(var, var2) 
Z <- c(z, z2) 
Significant <- c(sig, sig2) 
Pattern <- c(pattern, pattern2)


# make table dataframe holding column variables
data.for.table = data.frame(variables, moransI, expectedI, variance, Z, Significant, Pattern)


# create a graphical object using the table dataframe
# set the caption and table design specifics
table <- tableGrob(data.for.table, rows = c("","")) 
tCaption <- textGrob("Table 1: Global Moran's I statistics for Median Income and Refugee Population", 
                     just = "centre", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table <- gtable_add_rows(table, 
                         heights = grobHeight(t1Caption) + padding, 
                         pos = 0)

table <- gtable_add_grob(table,
                         tCaption, t = 1, l = 2, r = ncol(data.for.table1) + 1)


# generate the table
grid.arrange(table, newpage = TRUE)







#Printing the Tables (You can use the same setup for printing other types of objects (see ?png))
png("Output_Table1.png") #Create an object to print the table to
grid.arrange(table1, newpage = TRUE)
dev.off() #Print table










  
######################## 
#### LOCAL MORANS I ####
########################

# run a local Moran's I test
lisa.test <- localmoran(crd.data$MdInc, crd.lw)

# store results in new columns in the dataset
crd.data$Ii <- lisa.test[,1]
crd.data$E.Ii<- lisa.test[,2]
crd.data$Var.Ii<- lisa.test[,3]
crd.data$Z.Ii<- lisa.test[,4]
crd.data$P<- lisa.test[,5]


# map the local Moran's i results
# determine Moran's I legend parameters
minIi <- min(crd.data$Ii, na.rm=T)
maxIi <- max(crd.data$Ii, na.rm=T)
expIi <- lisa.test[1,2]

# create map
map_LISA <- tm_shape(crd.data) + 
  tm_polygons(col = "Ii", 
              title = "Local Moran's I", 
              style = "fixed",
              breaks = c(minIi, expIi-0.005, expIi+0.005, maxIi),
              labels = c("Negative Autocorrelation", "Random",
                         "Positive Autocorrelation"),
              midpoint = NA,
              palette = "BrBG", n = 3) +
  tm_layout(title = "Local Moran's I for Median Income in Victoria")
map_LISA


# map the p-values to visualize significance
# correct p-values for a two-tail significance test
crd.data$P2Tail <- crd.data$P * 2

# determine p-value legend parameters
minP <- min(crd.data$P2Tail, na.rm=T)
maxP <- max(crd.data$P2Tail, na.rm=T)

#create map
map_PVals <- tm_shape(crd.data) + 
  tm_polygons(col = "P2Tail", 
              title = "P Values", 
              style = "fixed",
              breaks = c(minP, 0.05, maxP),
              labels = c("p < 0.05", "p > 0.05"),
              palette = "RdGy", n = 3) +
  tm_layout(title = "P values from Median Income Local Moran's I Test")
map_PVals


# graph Local Moran's I results in a scatterplot
plot.MdInc <- moran.plot(crd.data$MdInc, crd.lw, 
        zero.policy=NULL, spChk=NULL, labels=NULL, 
        main= "Scatterplot of Median Income Autocorrelation by Census Tracts", 
        xlab="Median Income", ylab="Spatially Lagged Median Income", quiet=NULL)









######################
# REFUGEE LOCAL MORAN
####################


lisa.test2 <- localmoran(crd.data$RefIm, crd.lw)

crd.data$Ii2 <- lisa.test2[,1]
crd.data$E.Ii2<- lisa.test2[,2]
crd.data$Var.Ii2<- lisa.test2[,3]
crd.data$Z.Ii2<- lisa.test2[,4]
crd.data$P2<- lisa.test2[,5]


# map the local Moran's i results
# determine Moran's I legend parameters
minIi2 <- min(crd.data$Ii2, na.rm=T)
maxIi2 <- max(crd.data$Ii2, na.rm=T)
expIi2 <- lisa.test2[1,2]

# create map
map_LISA2 <- tm_shape(crd.data) + 
  tm_polygons(col = "Ii2", 
              title = "Local Moran's I", 
              style = "fixed",
              breaks = c(minIi2, expIi2-0.005, expIi2+0.005, maxIi2),
              labels = c("Negative Autocorrelation", "Random",
                         "Positive Autocorrelation"),
              midpoint = NA,
              palette = "BrBG", n = 3) +
  tm_layout(title = "Local Moran's I for Refugee Population in Victoria")
map_LISA2


# map the p-values to visualize significance
# correct p-values for a two-tail significance test
crd.data$P22Tail <- crd.data$P2 * 2

# determine p-value legend parameters
minP2 <- min(crd.data$P22Tail, na.rm=T)
maxP2 <- max(crd.data$P22Tail, na.rm=T)

map_PVals2 <- tm_shape(crd.data) + 
  tm_polygons(col = "P22Tail", 
              title = "P Values", 
              style = "fixed",
              breaks = c(minP2, 0.05, maxP2),
              labels = c("p < 0.05", "p > 0.05"),
              palette = "RdGy", n = 2) +
  tm_layout(title = "P values from Refugee Population Local Moran's I Test")
map_PVals2





###############################
### MAPPING LOCAL MORAN'S I ###
###############################

#mapping the local moran's i using jenks natural breaks classification again

tmaptools::palette_explorer() #Tool for selecting pallettes

map_LISA <- tm_shape(crd.data) + 
  tm_polygons(col = "Ii", 
              title = "Local Moran's I", 
              style = "jenks", 
              palette = "BrBG", n = 3) +
  tm_layout(title = "Local Moran's I for Median Income in Victoria")


map_LISA



map_LISA2 <- tm_shape(crd.data) + 
  tm_polygons(col = "Ii2", 
              title = "Local Moran's I", 
              style = "jenks", 
              palette = "BrBG", n = 5) +
  tm_layout(title = "Local Moran's I for Refugee Population in Victoria")


map_LISA2


#########################
### GRAPH THE RESULTS ###
#########################

# make a scatterplot giving median income values and the weighted matrix
# positive slope indicates positive autocorrelation.
# degree of slope indicates the strnegth of the correlation
# the little numbered targets are significant results having a pull on the line of best fit

plot.MdInc <- moran.plot(crd.data$MdInc, crd.lw, zero.policy=NULL, spChk=NULL, labels=NULL, 
           main= "Scatterplot of Median Income Autocorrelation by Census Tracts", xlab="Median Income", 
           ylab="Spatially Lagged Median Income", quiet=NULL)


plot.Ref <- moran.plot(crd.data$RefIm, crd.lw, zero.policy=NULL, spChk=NULL, labels=NULL, 
            xlab="Refugee Population", ylab="Spatially Lagged Refugee Population", 
            quiet=NULL, main= "Scatterplot of Refugee Population Autocorrelation in Victoria CRD")


########################





