#!/usr/bin/Rscript
###########################################################################
### [V] Classification supervisee base sur les pixeles en images 
### drone puis exportation des donnees sur des rayons de 25 a 400 m
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
###########################################################################

###########################################################################
### [0] CHARGEMENT DES PACKAGES ET REPERTOIRE DE TRAVAIL
###########################################################################

require("sp")
### Pebesma, E.J., R.S. Bivand, 2005. Classes and methods for spatial
### data in R. R News 5 (2), https://cran.r-project.org/doc/Rnews/.
### Roger S. Bivand, Edzer Pebesma, Virgilio Gomez-Rubio, 2013. Applied
### spatial data analysis with R, Second edition. Springer, NY.
### http://www.asdar-book.org/
require("RStoolbox")
### Benjamin Leutner and Ned Horning (2017). RStoolbox: Tools for Remote
### Sensing Data Analysis. R package version 0.1.10.
### https://CRAN.R-project.org/package=RStoolbox
require("raster")
### Robert J. Hijmans (2016). raster: Geographic Data Analysis and
### Modeling. R package version 2.5-8.
### https://CRAN.R-project.org/package=raster
require("randomForest")
### A. Liaw and M. Wiener (2002). Classification and Regression by
### randomForest. R News 2(3), 18--22.
require("e1071")
### David Meyer, Evgenia Dimitriadou, Kurt Hornik, Andreas Weingessel and
### Friedrich Leisch (2017). e1071: Misc Functions of the Department of
### Statistics, Probability Theory Group (Formerly: E1071), TU Wien. R
### package version 1.6-8. https://CRAN.R-project.org/package=e1071

wd <- "/home/bioinfo/Documents/2018_M2_CamilaBenavides/SIG_QGIS/DRONE"
setwd(wd)
# list.files()

###########################################################################
### [1] CLASSIFICATION SUPERVISEE
###########################################################################

### 1.1. importer les rasters et les training shapefiles 
miRasters <- list(
	raster("P5_transparent_mosaic_group1.tif"), 
	raster("P7_transparent_mosaic_group1.tif"), 
	raster("P12_transparent_mosaic_group1.tif"),
	raster("P18_transparent_mosaic_group1.tif"), 
	raster("P38_transparent_mosaic_group1.tif"), 
	raster("P40_transparent_mosaic_group1.tif"),
	raster("P41_transparent_mosaic_group1.tif"), 
	raster("P42_transparent_mosaic_group1.tif")
)
miShps <- list( ### liste des ROI
	shapefile("ROI_P05_MC.shp"), 
	shapefile("ROI_P07_MC.shp"), 
	shapefile("ROI_P12_MC.shp"), 
	shapefile("ROI_P18_MC.shp"), 
	shapefile("ROI_P38_MC.shp"), 
	shapefile("ROI_P40_MC.shp"), 
	shapefile("ROI_P41_MC.shp"), 
	shapefile("ROI_P42_CM.shp")
)

# miRaster <- miRasters[[1]]
# miShapeFile <- miShps[[1]]

### 1.2. verifier visuellement la superposition des couches
# plot(miRaster)
# plot(miShapeFile, add = TRUE) 

### 1.3. Faire la classification supervisee et verifier visuellement
# miClassSup_MC <- superClass(
  # img = miRaster, 
  # trainData = miShapeFile, 
  # responseCol = 1) ####### pour une seule parcelle

# par(mfrow = c(1, 2))
# plot(miRaster)
# plot(miClassSup_MC$map)

listMiClassSup_MC <- lapply(1:length(miRasters), function(i){
	miClassSup_MC <- superClass(
		img = miRasters[[i]], 
		trainData = miShps[[i]], 
		responseCol = 1)
})

###########################################################################
### [2] EXTRACTION DES DONNES PAR RAYON
###########################################################################

pointsDF <- data.frame(
	parcelle = c("P5", "P7", "P12", "P18", "P38", "P40", "P41", "P42"), 
	longitude = c(-78.78509333, -78.80273167, -78.79859667, -78.63821333, 
		-78.588908, -78.743102, -78.730944, -78.76065), 
	latitude = c(-2.108118333, -2.088856667, -2.100291667, -0.963956667, 
		-0.907498, -0.824985, -0.815995, -2.008947)
)
x <- pointsDF[,2]
y <- pointsDF[,3]
d <- data.frame(lon = x, lat = y)
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+init=epsg:4326") # WGS 84
CRS_new <- CRS("+proj=utm +zone=17 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
d_ch1903 <- spTransform(d, CRS_new)

radiusR <- seq(from = 25, to = 400, by = 25)

RS_drone_pixSup_K <- lapply(1:length(listMiClassSup_MC), function(numClassif){
	myPoint <- d_ch1903[numClassif]
	indi <- sapply(radiusR, function(i){
		myRaster_k <- extract(listMiClassSup_MC[[numClassif]]$map,
		myPoint,
		buffer = i,
		fun = function(k){
			k1 <- length(k[k == 1]) / length(k)
			k2 <- length(k[k == 2]) / length(k)
			k3 <- length(k[k == 3]) / length(k)
			return(list(k1, k2, k3))
		},
		df = FALSE)
	})
	rownames(indi) <- paste0(as.character(pointsDF$parcelle[numClassif]), "_", 
		c("k1", "k2", "k3")) 
	colnames(indi) <- paste0("radius", radiusR)
	return(indi)
})

save(listMiClassSup_MC, file = "listMiClassSup_MC.RData")

RS_drone_pixSup_percentK_DF <- data.frame(do.call(rbind, RS_drone_pixSup_K))
rownames(RS_drone_pixSup_percentK_DF) <- paste0(
	rep(as.character(pointsDF$parcelle), each = 3), 
	"_", 
	c("k1", "k2", "k3") # ...kn:le nombre de categories que nous avons
)
  
save(RS_drone_pixSup_percentK_DF, file = "RS_drone_pixSup_percentK_DF.RData")
print(RS_drone_pixSup_percentK_DF)
write.table(RS_drone_pixSup_percentK_DF, "RS_drone_pixSup_percentK_DF.xls", 
            col = NA, sep = "\t", dec = ".") 
