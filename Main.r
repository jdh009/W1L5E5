### Jeroen Roelofs
### January 10 2015

## Load packages and functions.
library(downloader)
library(rgdal)
library(sp)
library(rgeos)

# Download data.
download("http://www.mapcruzin.com/download-shapefile/netherlands-places-shape.zip", "Data/netherlands-places-shape.zip", quiet = T, mode = "wb")
download("http://www.mapcruzin.com/download-shapefile/netherlands-railways-shape.zip", "Data/netherlands-railways-shape.zip", quiet = T, mode = "wb")

# Unpackage data.
unzip("Data/netherlands-places-shape.zip", exdir = "Data//NetherlandsPlaces") # could also use: Data <- unz(temp, filename)
unzip("Data/netherlands-railways-shape.zip", exdir = "Data/NetherlandsRailways") # could also use: Data <- unz(temp, filename)

# Load required data.
railways <- readOGR("Data/NetherlandsRailways/railways.shp", "railways")
places <- readOGR("Data/NetherlandsPlaces/places.shp", "places" )
industrial <- railways[railways$type == 'industrial',]

# Check if the loading worked
# plot(railways)
# plot(industrial)
# plot(places)

# define CRS object for RD projection.
transformRD <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889
                                      +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.2369,
                                      50.0087,465.658,-0.406857330322398,0.350732676542563,-1.8703473836068,
                                      4.0812 +units=m +no_defs")

# perform the coordinate transformation from WGS84 to RD.
placesRD <- spTransform(places, transformRD)
industrialRailwayRD <- spTransform(industrial, transformRD)

# Check if the transformation worked.
# plot(industrialRailwayRD)
# plot(placesRD)

#  Create 1000m buffer around industral railway. 
industrialrailbuffer <- gBuffer(industrialRailwayRD, byid=F, width=1000)

# check buffer.
# plot(buffer)
# plot(industrialRailwayRD, add=TRUE)

#Intersection buffer. 
i <- gIntersects(industrialrailbuffer, placesRD, byid = TRUE) # extract name
coordinates <- gIntersection(placesRD, industrialrailbuffer, byid = TRUE) #extract coordinates: 5973 buffer 136033 456446.2
place <- placesRD@data[i] #[1] "235861650" "Utrecht"   "city"      "100000" 
data <- data.frame(Name = place[2], Population = place[4])
# placename <- place[2] #[1] "Utrecht"
utrecht <- SpatialPointsDataFrame(coords = coordinates, data, proj4string=transformRD, match.ID = FALSE)

# Assigning Coordinates.
# x <- coordinates@coords[1]
# x <- as.integer(x)
# y <- coordinates@coords[2]
# y <- as.integer(y)
# coordinates  <- coordinates(c(x,y))

# print all
plot(industrialrailbuffer, axes = TRUE, col='peachpuff2')
plot(industrialRailwayRD, add = TRUE)
plot(utrecht, cex = 1.5, add = TRUE)
mtext(side=2, "Latitude", line=2.5)
mtext(side=1, "Longitude", line=2.5)
text(utrecht@coords, labels=as.character(utrecht$Name), cex=1.1, font=7, pos=2, offset=0.8)
# invisible(text(getSpPPolygonsLabptSlots(placename), labels = place.label[2], cex = 1.2, col = "white", font = 9))

# Point 5: write down the name of the city and the population of that city as one comment at the end of the script.
print(paste("City :", place[2])) #Utrecht
print(paste("Population :", place[4])) #100.000