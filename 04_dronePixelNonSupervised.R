#!/usr/bin/Rscript
###############################################################################
### [IV] Classification non supervise basee sur les pixels et images drone
### puis exportation des donnees sur des rayons de 50 a 400 m
###
### IRD / Camila Benavides, Francois Rebaudo, 2018
### License: Creative Commons CC-BY-NC-SA
###
### To cite this script please use :
### Benavides C. and Rebaudo F. 2018; R script to analyze Remote Sensing 
### Imagery in landscape ecology; in: Étude du niveau de phytophagie sur la 
### culture de lupin (Lupinus mutabilis Sweet) dans les Andes équatoriennes 
### et sa relation à la configuration et composition du paysage, Master tesis, 
### AgroParisTech, France
###############################################################################

###############################################################################
### [0] CHARGEMENT DES PACKAGES ET REPERTOIRE DE TRAVAIL
###############################################################################

require("RStoolbox")
### Benjamin Leutner and Ned Horning (2017). RStoolbox: Tools for Remote
### Sensing Data Analysis. R package version 0.1.10.
### https://CRAN.R-project.org/package=RStoolbox
require("raster")
### Robert J. Hijmans (2016). raster: Geographic Data Analysis and
### Modeling. R package version 2.5-8.
### https://CRAN.R-project.org/package=raster
require("rgdal")
### Roger Bivand, Tim Keitt and Barry Rowlingson (2017). rgdal: Bindings
### for the 'Geospatial' Data Abstraction Library. R package version
### 1.2-15. https://CRAN.R-project.org/package=rgdal
require("dplyr")
### Hadley Wickham, Romain Francois, Lionel Henry and Kirill Müller
### (2017). dplyr: A Grammar of Data Manipulation. R package version
### 0.7.4. https://CRAN.R-project.org/package=dplyr
require("maptools")
### Roger Bivand and Nicholas Lewin-Koh (2017). maptools: Tools for
### Reading and Handling Spatial Objects. R package version 0.9-2.
### https://CRAN.R-project.org/package=maptools

wd <- "/home/bioinfo/Documents/2018_M2_CamilaBenavides/SIG_QGIS/DRONE/"
setwd(wd)
# list.files()

###############################################################################
### [1] CLASSIFICATION NON SUPERVISEE
###############################################################################

# temps de calcul d'environ une heure :
# le fichier .tif cest l'image a utiliser
# a faire une seule fois puis utiliser le fichier listRASTERUnClass5.RData
pixelUnClass_p5 <- unsuperClass(img = raster("P5_transparent_mosaic_group1.tif"), 
	nClasses = 3)
pixelUnClass_p7 <- unsuperClass(img = raster("P7_transparent_mosaic_group1.tif"), 
	nClasses = 3)
pixelUnClass_p12 <- unsuperClass(img = raster("P12_transparent_mosaic_group1.tif"), 
	nClasses = 3)
pixelUnClass_p18 <- unsuperClass(img = raster("P18_transparent_mosaic_group1.tif"), 
	nClasses = 3)
pixelUnClass_p38 <- unsuperClass(img = raster("P38_transparent_mosaic_group1.tif"), 
	nClasses = 3)
pixelUnClass_p40 <- unsuperClass(img = raster("P40_transparent_mosaic_group1.tif"), 
	nClasses = 3)
pixelUnClass_p41 <- unsuperClass(img = raster("P41_transparent_mosaic_group1.tif"), 
	nClasses = 3)
pixelUnClass_p42 <- unsuperClass(img = raster("P42_transparent_mosaic_group1.tif"), 
	nClasses = 3)

listRASTERUnClass5 <- list(pixelUnClass_p5, pixelUnClass_p7, pixelUnClass_p12, 
	pixelUnClass_p18, pixelUnClass_p38, pixelUnClass_p40, pixelUnClass_p41, 
	pixelUnClass_p42)

save(listRASTERUnClass5, file = "listRASTERUnClass5.RData")

rm("pixelUnClass_p5", "pixelUnClass_p7", "pixelUnClass_p12", 
     "pixelUnClass_p18", "pixelUnClass_p38", "pixelUnClass_p40", 
	 "pixelUnClass_p41", "pixelUnClass_p42")

load("listRASTERUnClass5.RData")

### representation graphique de la classification non supervisee
# colors <- rainbow(5)
# plotName <- paste0("P_", 1:8)
# par(mfrow = c(3, 3))
# for(i in 1:length(listRASTER)){
	# plot(listRASTER[[i]]$map, col = colors, 
		# legend = FALSE, 
		# axes = FALSE, 
		# box = FALSE, 
		# main = plotName[i])
# }

###############################################################################
### [2] EXTRACTION DES DONNES PAR RAYON (3 categories)
###############################################################################

################### valeurs de deviation (DS) ################################
pointsDF <- data.frame(
	parcelle = c("P5", "P7", "P12", "P18", "P38", "P40", "P41", "P42"), 
	longitude = c(-78.78509333, -78.80273167, -78.79859667, -78.63821333, 
		-78.588908, -78.743102, -78.730944, -78.76065), 
	latitude = c(-2.108118333, -2.088856667, -2.100291667, -0.963956667, 
		-0.907498, -0.824985, -0.815995, -2.008947)
)

radiusR <- seq(from = 25, to = 400, by = 25)

x <- pointsDF[,2]
y <- pointsDF[,3]
d <- data.frame(lon = x, lat = y)
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+init=epsg:4326") # WGS 84
CRS_new <- CRS("+proj=utm +zone=17 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
d_ch1903 <- spTransform(d, CRS_new)

### verifier visuellement la superposition des couches:
# plot(myRaster$map)
# points(d_ch1903)

### obtenir la DS aux differents rayons:
RS_drone_pixUnSup_sdK <- lapply(1:length(listRASTERUnClass5), function(myRaster){
	posPoints <- d_ch1903
	  myPoint <- posPoints[myRaster]
	  indi <- sapply(radiusR, function(i){
			myRaster_k <- extract(listRASTERUnClass5[[myRaster]]$map,
			myPoint,
			buffer = i,
			fun = function(k){
				k1 <- length(k[k == 1]) / length(k)
				k2 <- length(k[k == 2]) / length(k)
				k3 <- length(k[k == 3]) / length(k)
			  return(sd(c(k1, k2, k3), na.rm = TRUE))
			},
			df = FALSE)
		})
	rasname <- strsplit(as.character(listRASTERUnClass5[[myRaster]]$call)[2], 
		split = "\"")[[1]][2]
	names(indi) <- paste0(rasname, "_sd_", as.character(radiusR))
	return(indi)
})

RS_drone_pixUnSup_DF <- data.frame(do.call(cbind, RS_drone_pixUnSup_sdK))
rownames(RS_drone_pixUnSup_DF) <- paste0("radius", radiusR)
colnames(RS_drone_pixUnSup_DF) <- as.character(pointsDF$parcelle)
RS_drone_pixUnSup_DF_percentK <- RS_drone_pixUnSup_DF
save(RS_drone_pixUnSup_DF_percentK, file = "RS_drone_pixUnSup_DF_percentK.RData")

print(RS_drone_pixUnSup_DF_percentK)

############################ valeurs maximales ##################################

RS_drone_pixUnSup_maxK <- lapply(1:length(listRASTERUnClass5), function(myRaster){
  posPoints <- d_ch1903
  myPoint <- posPoints[myRaster]
  indi <- sapply(radiusR, function(i){
    myRaster_k <- extract(listRASTERUnClass5[[myRaster]]$map,
                          myPoint,
                          buffer = i,
                          fun = function(k){
                            k1 <- length(k[k == 1]) / length(k)
                            k2 <- length(k[k == 2]) / length(k)
                            k3 <- length(k[k == 3]) / length(k)
                            return(max(c(k1, k2, k3), na.rm = TRUE))
                          },
                          df = FALSE)
  })
  rasname <- strsplit(as.character(listRASTERUnClass5[[myRaster]]$call)[2], 
	split = "\"")[[1]][2]
  names(indi) <- paste0(rasname, "_sd_", as.character(radiusR))
  return(indi)
})

RS_drone_pixUnSup_DFmax <- data.frame(do.call(cbind, RS_drone_pixUnSup_maxK))
rownames(RS_drone_pixUnSup_DFmax) <- paste0("radius", radiusR)
colnames(RS_drone_pixUnSup_DFmax) <- as.character(pointsDF$parcelle)
RS_drone_pixUnSup_DF_maxK <- RS_drone_pixUnSup_DFmax

save(RS_drone_pixUnSup_DF_maxK, file = "RS_drone_pixUnSup_DF_maxK.RData")
print(RS_drone_pixUnSup_DF_maxK)

par(mfrow = c(1, 2))
boxplot(RS_drone_pixUnSup_DF_percentK, main = "sd")
boxplot(RS_drone_pixUnSup_DF_maxK, main = "max")
