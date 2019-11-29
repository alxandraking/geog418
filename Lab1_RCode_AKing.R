# Alexandra King
# V00827380
# Assignment 1 - Geog418 Geostatistics


#Install Libraries - don't have to do every time! (If you see "a package doesnt exist with that name"  then re-install)
install.packages("rgdal")
install.packages("lubridate")
install.packages("tidyverse")
install.packages("grid")
install.packages("gridExtra")
install.packages("gtable")
install.packages("maps")
install.packages("bcmaps")
install.packages("tmap")
install.packages("fBasics")
install.packages('bcmapsdata', repos='https://bcgov.github.io/drat/')


#Load Libraries - if error "function not found" pops up, your libraries aren't loaded
library("rgdal")
library("lubridate")
library("tidyverse")
library("grid")
library("gridExtra")
library("gtable")
library("maps")
library("bcmaps")
library("tmap")
library("fBasics")
library('bcmapsdata', repos='https://bcgov.github.io/drat/')

#####
#Set working directory
dir <- "Z:\\Geog418\\Assignment1"
setwd(dir)
getwd()

#####
#Read in data and extract attribute table
shp <- readOGR(dsn = ".\\Working", layer = "prot_current_fire_points") #read in shp file from current (".") working directory

#here we create a variable called 'df' which is the attribute table stored in a dataframe
df <- shp@data #Extract the attribute table as a dataframe, store in variable
class(df) #ensure new variable is a dataframe
attach(df) #attach dataset

View(df) #allows to view attribute table

#####
#Inspect data
names(df) #see column names
head(df) #see first 6 rows of data

range(df$FIRE_YEAR) #How many years of data is there? dollar sign means you're calling on that datafram for that column

nrow(subset(df, FIRE_YEAR == 2018)) #How many observations in 2018?
nrow(subset(df, FIRE_YEAR == 2019)) #How many observations in 2019?

#####
#Subset 2019 data and calculate descriptive statistics
df_2019 <- df[which(df$FIRE_YEAR == 2019),] #subset only 2019 data

#Mean
meanPop <- mean(df_2019$CURRENT_SI, na.rm = TRUE) #Use na.rm = TRUE to ignore NA values in calculation
meanPop <- round(meanPop, digits=3)

df_2019$IGNITION_D <- as.Date(df_2019$IGNITION_D) #Convert from factor to date string
df_2019$IGN_DJUL <- yday(df_2019$IGNITION_D) #Make new column with Julian day for the ignition date
meanSummer <- mean(subset(df_2019, IGN_DJUL >= 182 & IGN_DJUL <= 243)$CURRENT_SI) #Calculate the mean fire size between July 1 (182) and Aug 31 (243)
meanSummer <- round(meanSummer, digits = 3)

#Standard Deviation
sdPop <- sd(df_2019$CURRENT_SI, na.rm = TRUE) #Calculate the SD, ignoring NA values
sdPop <- round(sdPop, digits = 3)
sdSummer <- sd(subset(df_2019, IGN_DJUL >= 182 & IGN_DJUL <= 243)$CURRENT_SI) #Calculate the SD, ignoring NA values only for the summer months
sdSummer <- round(sdSummer, digits = 3)

#Mode
modePop <- as.numeric(names(sort(table(df_2019$CURRENT_SI), decreasing = TRUE))[1]) #make frequency table of fire size variable and sort it in desending order and extract the first row (Most Frequent)
modePop <- round(modePop, digits = 3)

#It will be cleaner if we use a new variable for the summer data
df_Summer <- subset(df_2019, IGN_DJUL >= 182 & IGN_DJUL <= 243) #Make new variable to hold summer data
#make frequency table of fire size variable and sort it in desending order and extract the first row (Most Frequent)
modeSummer <- as.numeric(names(sort(table(df_Summer$CURRENT_SI), decreasing = TRUE))[1])
modeSummer <- round(modeSummer, digits = 3)

#Median
medPop <- median(df_2019$CURRENT_SI, na.rm = TRUE)
medPop <- round(medPop, digits = 3)
medSummer <- median(df_Summer$CURRENT_SI, na.rm = TRUE)
medSummer <- round(medSummer, digits = 3)
  
#Skewness
skewPop <- skewness(df_2019$CURRENT_SI, na.rm = TRUE)[1]
skewPop <- round(skewPop, digits = 3)
skewSummer <- skewness(df_Summer$CURRENT_SI, na.rm = TRUE)[1]
skewSummer <- round(skewSummer, digits = 3)
  
#Kurtosis
kurtPop <- kurtosis(df_2019$CURRENT_SI, na.rm = TRUE)[1]
kurtPop <- round(kurtPop, digits = 3)
kurtSummer <- kurtosis(df_Summer$CURRENT_SI, na.rm = TRUE)[1]
kurtSummer <- round(kurtSummer, digits = 3)

#CoV - a measure of you dispersion compared to your mean. its the SD divided by Mean x 100. a RELATIVE measurement to your mean
CoVPop <- (sdPop / meanPop) * 100
CoVPop <- round(CoVPop, digits = 3)
CoVSummer <- (sdSummer / meanSummer) * 100
CoVSummer <- round(CoVSummer, digits = 3)

#Normal distribution test
#Shapiro test = a test to see if your distribution is normal or not
normPop_PVAL <- shapiro.test(df_2019$CURRENT_SI)$p.value
normPop_PVAL <- round(normPop_PVAL, digits = 3)
normSummer_PVAL <- shapiro.test(df_Summer$CURRENT_SI)$p.value
normSummer_PVAL <- round(normSummer_PVAL, digits = 3)


#####
#Create a table of descriptive stats

samples = c("Population", "Summer") #Create an object for the labels
means = c(meanPop, meanSummer) #Create an object for the means
standardDeviation = c(sdPop, sdSummer) #Create an object for the standard deviations
medians = c(medPop, medSummer) #Create an object for the medians
modes <- c(modePop, modeSummer) #Create an object for the modes
skewness <- c(skewPop, skewSummer) #Create an object for the skewness
kurtosis <- c(kurtPop, kurtSummer) #Create an object for the kurtosis
CoV <- c(CoVPop, CoVSummer) #Create an object for the CoV
normality <- c(normPop_PVAL, normSummer_PVAL) #Create an object for the normality PVALUE

data.for.table1 = data.frame(samples, means, standardDeviation, medians, modes)
data.for.table2 = data.frame(samples, skewness, kurtosis, CoV, normality)


#Make table 1
table1 <- tableGrob(data.for.table1, rows = c("","")) #make a table "Graphical Object" (GrOb) 
t1Caption <- textGrob("Table 1: Descriptive statistics of mean, standard dev., median, and 
                    mode for 2019 wildfire data and summer months of July and August.", just = 0.54, gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table1 <- gtable_add_rows(table1, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0)

table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.table1) + 1) #because 0 counts as 1 in a index.

#Make table 2
table2 <- tableGrob(data.for.table2, rows = c("",""))
t2Caption <- textGrob("Table 2: Descriptive statistics of skewness, kurtosis, C. of V., and 
                    Shapiro test for 2019 wildfire data and summer months of July & Aug.", just = 0.57, gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table2 <- gtable_add_rows(table2, 
                          heights = grobHeight(t2Caption) + padding, 
                          pos = 0)

table2 <- gtable_add_grob(table2,
                          t2Caption, t = 1, l = 2, r = ncol(data.for.table2) + 1)


#Generate the Tables
grid.arrange(table1, newpage = TRUE)
grid.arrange(table2, newpage = TRUE)

#Printing the Tables (You can use the same setup for printing other types of objects (see ?png))
png("Output_Table1.png") #Create an object to print the table to
grid.arrange(table1, newpage = TRUE)
dev.off() #Print table

png("Output_Table2.png") #Create an object to print the table to
grid.arrange(table2, newpage = TRUE) #Create table
dev.off()

#####
#Create and Print a histogram
#png("Output_Histogram.png")
#hist(df_2019$CURRENT_SI, breaks = 30, main = "TITLE OF YOUR HISTOGRAM", xlab = "LABEL FOR AXIS") #Base R style
#dev.off()

#Creating Histogram 1: all 2019 wildfire data

histogram <- ggplot(df_2019, aes(x = CURRENT_SI)) + #Create new GGplot object with data attached and fire size mapped to X axis
  geom_histogram(bins = 30, color = "black", fill = "white") + #make histogram with 30 bins, black outline, white fill
  labs(title = "Frequency of Fire Sizes for 2019", x = "Size of Wildfire Incidents (hectares)", y = "Frequency", caption = "Figure 1: Histogram showing the frequency of wildfire sizes for all 2019 wildfire data.") + #label plot, x axis, y axis
  theme_classic() + #set the theme to classic (removes background and borders etc.)
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.caption = element_text(hjust = 0.5)) + #set title to center and bold
  scale_y_continuous(breaks = seq(0, 700, by = 100)) # set y axis labels to 0 - 700 incrimenting by 100


#geom - giving ggplot the geometry of a histogram. bins - 30. color - black outline with white fill. next line is labels. 
#last line  setting scale of y-axis. useful for making multiple figures, all axis being same scale makes comparisons easier.

#Creating Histogram 2: summer months of July and August

histogram2 <- ggplot(df_Summer, aes(x = CURRENT_SI)) + #Create new GGplot object with data attached and fire size mapped to X axis
  geom_histogram(bins = 30, color = "black", fill = "white") + #make histogram with 30 bins, black outline, white fill
  labs(title = "Frequency of Fire Sizes for July & August 2019", x = "Size of Wildfire Incidents (hectares)", y = "Frequency", caption = "Figure 2: Histogram showing the frequency of wildfire sizes for the summer months of July and August.") + #label plot, x axis, y axis
  theme_classic() + #set the theme to classic (removes background and borders etc.)
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.caption = element_text(hjust = 0.5)) + #set title to center and bold
  scale_y_continuous(breaks = seq(0, 700, by = 100)) # set y axis labels to 0 - 700 incrimenting by 100

#Printing Histogram 1
png("Output_Histogram_2019.png")
histogram
dev.off()

#Printing Histogram 2
png("Output_Histogram_Summer.png")
histogram2
dev.off()



#####
#Creating bar graph
df_2019$IGN_Month <- month(df_2019$IGNITION_D, label = TRUE, abbr = TRUE) #create new column with month of ignition


sumMar = sum(subset(df_2019, IGN_Month == "Mar")$CURRENT_SI, na.rm = TRUE) #create new object for March
sumApr = sum(subset(df_2019, IGN_Month == "Apr")$CURRENT_SI, na.rm = TRUE) #create new object for April
sumMay = sum(subset(df_2019, IGN_Month == "May")$CURRENT_SI, na.rm = TRUE) #create new object for May
sumJun = sum(subset(df_2019, IGN_Month == "Jun")$CURRENT_SI, na.rm = TRUE) #create new object for June
sumJul = sum(subset(df_2019, IGN_Month == "Jul")$CURRENT_SI, na.rm = TRUE) #create new object for July
sumAug = sum(subset(df_2019, IGN_Month == "Aug")$CURRENT_SI, na.rm = TRUE) #create new object for August
sumSep = sum(subset(df_2019, IGN_Month == "Sep")$CURRENT_SI, na.rm = TRUE) #create new object for September
months = c("Mar","Apr","May","Jun","Jul", "Aug", "Sep")  #Create labels for the bar graph

png("Output_BarGraph.png") #Create an object to print the bar graph 
barplot(c(sumMar,sumApr,sumMay, sumJun, sumJul, sumAug, sumSep), names.arg = months, 
        main = "Total Hectares Burned by Month in 2019", ylab = "Hectares Burned", xlab = "Months") #Create the bar graph
dev.off() #Print bar graph

#Total Size by Month GGPLOT
barGraph <- df_2019 %>% #store graph in bargraph variable and pass data frame as first argument in next line
  group_by(IGN_Month) %>% #use data frame and group by month and pass to first argument in next line
  summarise(sumSize = sum(CURRENT_SI, na.rm = TRUE)) %>% #sum up the total fire size for each month and pass to GGplot
  ggplot(aes(x = IGN_Month, y = sumSize)) + #make new GGPLOT with summary as data and month and total fire size as x and y
  geom_bar(stat = "identity") + #make bar chart with the Y values from the data (identity)
  labs(title = "Total Hectares Burned by Month in 2019", x = "Months", y = "Hectares Burned", caption = "Figure 2: Bar graph showing the total hectares burned in each month of the 2019 fire season.") + #label plot, x axis, y axis
  theme_classic() + #set the theme to classic (removes background and borders etc.)
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.caption = element_text(hjust = 0.5)) #set title to center and bold
barGraph

png("Output_BarGraph_GG.png")
barGraph
dev.off()


#####
#Lets put it all together
pdf("Lab_1_Figures_and_Tables.pdf", onefile = TRUE)
grid.arrange(table1, newpage = TRUE)
grid.arrange(table2, newpage = TRUE)
histogram
barGraph
dev.off()



#####
#Creating maps
#First example
# bc <- as_Spatial(bc_neighbours()) #Get shp of BC bounds
# bc <- spTransform(bc, CRS("+init=epsg:4326")) #project to WGS84 geographic (Lat/Long)
# 
# bc <- bc[which(bc$name == "British Columbia" ),] #Extract just the BC province
# 
# png("FirstMap.png")
#   map(bc, fill = TRUE, col = "white", bg = "lightblue", ylim = c(40, 70)) #make map of province
#   points(df_2019$LONGITUDE ,df_2019$LATITUDE , col = "red", pch = 16) #add fire points
# dev.off()

#####
#Making Maps with tm package
#Make spatial object (Spatial points dataframe) out of data
coords <- df_2019[, c("LONGITUDE", "LATITUDE")] #Store coordinates in new object
crs <- CRS("+init=epsg:4326") #store the coordinate system (CRS) in a new object

firePoints <- SpatialPointsDataFrame(coords = coords, data = df_2019, proj4string = crs) #Make new spatial Points object using coodinates, data, and projection

#Calculating the Mean Centre coordinates by using the Latitude and Longitude columns

meanLat <- mean(df_2019$LATITUDE) 
meanLong <- mean(df_2019$LONGITUDE) 

meanCentre <- data.frame(name = "Mean Centre", long = meanLong, lat = meanLat)
meanCoords <- meanCentre[,c("long", "lat")]
meanCRS <- CRS("+init=epsg:4326")
meanPoint <- SpatialPointsDataFrame(coords = meanCoords, data = meanCentre, proj4string = crs)

map_TM <- tm_shape(bc) + #make the main shape
  tm_fill(col = "gray50") +  #fill polygons
  tm_shape(firePoints) +
  tm_symbols(col = "red", alpha = 0.3) +
  tm_shape(meanPoint) +
  tm_symbols(col = "yellow", alpha = 0.8) +
  tm_add_legend(type = "symbol", labels = c("Fires", "MeanCentre"), col = c(adjustcolor("red", alpha.f = 0.3),adjustcolor("yellow", alpha.f = 0.3)), shape = c(19,19)) +
  tm_layout(title = "BC Fire Locations 2019", title.position = c("LEFT", "BOTTOM"), legend.position = c("RIGHT", "TOP"))

map_TM



png("TmMap.png")
  map_TM
dev.off()
