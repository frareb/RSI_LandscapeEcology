#!/usr/bin/Rscript
###########################################################################
### [VI] DSM DTM 
### surface puis exportation des donnees sur des rayons de 25 a 400 m
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

require("raster")
### Robert J. Hijmans (2016). raster: Geographic Data Analysis and
### Modeling. R package version 2.5-8.
### https://CRAN.R-project.org/package=raster

wd <- "/home/bioinfo/Documents/2018_M2_CamilaBenavides/SIG_QGIS/"
setwd(wd)
# list.files()

###########################################################################
### [1] Terrain characteristics from package raster
### from ?terrain:
### Compute slope, aspect and other terrain characteristics from a raster 
### with elevation data. The elevation data should be in map units (typically 
### meter) for projected (planar) raster data. They should be in meters when 
### the coordinate reference system (CRS) is longitude/latitude.
###
### opt: Character vector containing one or more of these options: slope, 
### aspect, TPI, TRI, roughness, flowdir (see Details)
###########################################################################

wd <- "/home/bioinfo/Documents/2018_M2_CamilaBenavides/SIG_QGIS/DRONE"
setwd(wd)
list.files()

###########################################################################
### [1] CHARGEMENT DES DSM
###########################################################################
myRasterP5_DSM <- raster("P5_dsm.tif")
myRasterP7_DSM <- raster("P7_dsm.tif")
myRasterP12_DSM <- raster("P12_dsm.tif")
myRasterP18_DSM <- raster("P18_dsm.tif")
myRasterP38_DSM <- raster("P38_dsm.tif")
myRasterP40_DSM <- raster("P40_dsm.tif")
myRasterP41_DSM <- raster("P41_dsm.tif")
myRasterP42_DSM <- raster("P42_dsm.tif")

listRaster <- list(myRasterP5_DSM, myRasterP7_DSM, myRasterP12_DSM, 
	myRasterP18_DSM, myRasterP38_DSM, myRasterP40_DSM, myRasterP41_DSM, 
	myRasterP42_DSM) 

###########################################################################
### [2] CALCUL DES INDICES
###########################################################################

listTerrainSlope <- lapply(listRaster, function(i){terrain(x = i, opt = "slope")})
listTerrainAspect <- lapply(listRaster, function(i){terrain(x = i, opt = "aspect")})
listTerrainRoughness <- lapply(listRaster, function(i){terrain(x = i, opt = "roughness")})
#listTerrainFlowdir <- lapply(listRaster, function(i){terrain(x = i, opt = "flowdir")})
# flowdir: pas relevant pour notre etude
listTerrainHillShade <- lapply(1:8, function(i){
	hillShade(slope = listTerrainSlope[[i]], aspect = listTerrainAspect[[i]])
})

### exemple de visualisation sur la parcelle 5
# mySlope <- terrain(x = myRasterP5_DSM, opt = "slope")
# myAspect <- terrain(x = myRasterP5_DSM, opt = "aspect")
# myTPI <- terrain(x = myRasterP5_DSM, opt = "TPI") # Topographic Position Index
# myTRI <- terrain(x = myRasterP5_DSM, opt = "TRI") # Terrain Ruggedness Index
# myRoughness <- terrain(x = myRasterP5_DSM, opt = "roughness")
# myFlowdir <- terrain(x = myRasterP5_DSM, opt = "flowdir")
# myHillShade <- hillShade(slope = mySlope, aspect = myAspect)

# pdf(file = "plot_terrainIndices.pdf")
	# plot(mySlope, main = "Slope")
	# plot(myAspect, main = "Aspect")
	# # plot(myTPI, main = "TPI")
	# # plot(myTRI, main = "TRI")
	# plot(myRoughness, main = "Roughness", zlim = c(0, 0.15))
	# plot(myFlowdir, main = "Flowdir")
	# plot(myHillShade, main = "HillShade")
# dev.off()

###############################################################################
### [3] EXTRACTION DES DONNES PAR RAYON
###############################################################################

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

### SLOPE
RS_drone_DSM_meanSlope <- lapply(1:length(listTerrainSlope), function(myRaster){
	posPoints <- d_ch1903
	myPoint <- posPoints[myRaster]
	indi <- sapply(radiusR, function(i){
	  myRaster_k <- extract(listTerrainSlope[[myRaster]],
		myPoint,
		buffer = i,
		fun = mean,
		df = FALSE)
	})
	rasname <- names(listTerrainSlope[[myRaster]])
	names(indi) <- paste0(rasname, "_mean_", as.character(radiusR))
	return(indi)
})

RS_drone_DSM_meanSlope_DF <- data.frame(do.call(cbind, RS_drone_DSM_meanSlope))
rownames(RS_drone_DSM_meanSlope_DF) <- paste0("radius", radiusR)
colnames(RS_drone_DSM_meanSlope_DF) <- as.character(pointsDF$parcelle)
save(RS_drone_DSM_meanSlope_DF, file = "RS_drone_DSM_meanSlope_DF.RData")
print(RS_drone_DSM_meanSlope_DF)

### fonction generalisee
getRSDroneDSM_mean <- function(listTerrainX){
	RS_drone_DSM_meanXXX <- lapply(1:length(listTerrainX), function(myRaster){
		posPoints <- d_ch1903
		  myPoint <- posPoints[myRaster]
		  indi <- sapply(radiusR, function(i){
				myRaster_k <- extract(listTerrainX[[myRaster]],
				myPoint,
				buffer = i,
				fun = mean,
				df = FALSE)
			})
		rasname <- names(listTerrainX[[myRaster]])
		names(indi) <- paste0(rasname, "_mean_", as.character(radiusR))
		return(indi)
	})
	RS_drone_DSM_meanXXX_DF <- data.frame(do.call(cbind, RS_drone_DSM_meanXXX))
	rownames(RS_drone_DSM_meanXXX_DF) <- paste0("radius", radiusR)
	colnames(RS_drone_DSM_meanXXX_DF) <- as.character(pointsDF$parcelle)
	return(RS_drone_DSM_meanXXX_DF)
}

### ASPECT
RS_drone_DSM_meanAspect_DF <- getRSDroneDSM_mean(listTerrainX = listTerrainAspect)
save(RS_drone_DSM_meanAspect_DF, file = "RS_drone_DSM_meanAspect_DF.RData")
print(RS_drone_DSM_meanAspect_DF)

### ROUGHNESS
RS_drone_DSM_meanRoughness_DF <- getRSDroneDSM_mean(listTerrainX = listTerrainRoughness)
save(RS_drone_DSM_meanRoughness_DF, file = "RS_drone_DSM_meanRoughness_DF.RData")
print(RS_drone_DSM_meanRoughness_DF)

### FLOWDIR
RS_drone_DSM_meanFlowdir_DF <- getRSDroneDSM_mean(listTerrainX = listTerrainFlowdir)
save(RS_drone_DSM_meanFlowdir_DF, file = "RS_drone_DSM_meanFlowdir_DF.RData")
print(RS_drone_DSM_meanFlowdir_DF)

### HILLSHADE
RS_drone_DSM_meanHillShade_DF <- getRSDroneDSM_mean(listTerrainX = listTerrainHillShade)
save(RS_drone_DSM_meanHillShade_DF, file = "RS_drone_DSM_meanHillShade_DF.RData")
print(RS_drone_DSM_meanHillShade_DF)
