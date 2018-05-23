#!/usr/bin/Rscript
###############################################################################
### [II] Script pour le calcul des indices dans un rayon de 25 a 400 m
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

wd <- "/home/bioinfo/Documents/2018_M2_CamilaBenavides/SIG_QGIS/"
setwd(wd)
# list.files()

###############################################################################
### [1] CALCUL DES INDICES AUTOUR DES ZONES D INTERET
###############################################################################

### [1.0] Charger les coordonnees des points d interet
###############################################################################
pointsDF <- data.frame(
	parcelle = c("P5", "P7", "P12", "P18", "P38", "P40", "P41", "P42"), 
	longitude = c(-78.78509333, -78.80273167, -78.79859667, -78.63821333, 
		-78.588908, -78.743102, -78.730944, -78.76065), 
	latitude = c(-2.108118333, -2.088856667, -2.100291667, -0.963956667, 
		-0.907498, -0.824985, -0.815995, -2.008947)
)

### [1.1] Charger les raster des indices au niveau du paysage
###############################################################################
# setwd("..")
tassCap_brigh <- raster("TassCap_brightness.tif")
tassCap_green <- raster("TassCap_greenness.tif")
tassCap_wetne <- raster("TassCap_wetness.tif")
specInd_CTVI <- raster("specInd_CTVI.tif")
specInd_DVI <- raster("specInd_DVI.tif")
specInd_GEMI <- raster("specInd_GEMI.tif")
specInd_GNDVI <- raster("specInd_GNDVI.tif")
specInd_MNDWI <- raster("specInd_MNDWI.tif")
specInd_MSAVI <- raster("specInd_MSAVI.tif")
specInd_MSAVI2 <- raster("specInd_MSAVI2.tif")
specInd_NBRI <- raster("specInd_NBRI.tif")
specInd_NDVI <- raster("specInd_NDVI.tif")
specInd_NDWI <- raster("specInd_NDWI.tif")
specInd_NDWI2 <- raster("specInd_NDWI2.tif")
specInd_NRVI <- raster("specInd_NRVI.tif")
specInd_RVI <- raster("specInd_RVI.tif")
specInd_SATVI <- raster("specInd_SATVI.tif")
specInd_SAVI <- raster("specInd_SAVI.tif")
specInd_SLAVI <- raster("specInd_SLAVI.tif")
specInd_SR <- raster("specInd_SR.tif")
specInd_TVI <- raster("specInd_TVI.tif")
specInd_TTVI <- raster("specInd_TTVI.tif")
specInd_WDVI <- raster("specInd_WDVI.tif")

listRASTER <- list(tassCap_brigh, tassCap_green, tassCap_wetne, 
	specInd_CTVI, specInd_DVI, specInd_GEMI, specInd_GNDVI, 
	specInd_MNDWI, specInd_MSAVI, specInd_MSAVI2, specInd_NBRI, 
	specInd_NDVI, specInd_NDWI, specInd_NDWI2, specInd_NRVI, 
	specInd_RVI, specInd_SATVI, specInd_SAVI, specInd_SLAVI, 
	specInd_SR, specInd_TVI, specInd_TTVI, specInd_WDVI)

### [1.2] Calcul des indices dans une zone tampon circulaire
###############################################################################
radiusR <- seq(from = 25, to = 400, by = 25)

x <- pointsDF[,2]
y <- pointsDF[,3]
d <- data.frame(lon = x, lat = y)
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+init=epsg:4326") # WGS 84
CRS_new <- CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
d_ch1903 <- spTransform(d, CRS_new)

RS_satellite <- lapply(listRASTER, function(myRaster){
	posPointsName <- as.character(pointsDF[,1])
	posPoints <- d_ch1903
	indi <- sapply(radiusR, function(i){
		myRaster_mean <- extract(myRaster,
			posPoints,
			buffer = i,
			fun = mean,
			df = FALSE)
	})
	colnames(indi) <- paste0(names(myRaster), "_mean_", as.character(radiusR))
	rownames(indi) <- posPointsName
	return(indi)
})
RS_satellite_mean_DF <- data.frame(do.call(cbind, RS_satellite))
RS_satellite_mean_DF$ID <- rownames(RS_satellite_mean_DF)
## sauvegarde des resultats dans des fichiers
save(RS_satellite_mean_DF, file = "RS_satellite_mean_DF.RData")
write.table(t(RS_satellite_mean_DF), "RS_satellite_mean_DF.xls", 
	col = NA, sep = "\t", dec = ",") 

