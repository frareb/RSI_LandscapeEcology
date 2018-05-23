#!/usr/bin/Rscript
################################################################################
### [IX] Script des analyses: indices Drone configuration vs. dommage
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
################################################################################

################################################################################
### [0] CHARGEMENT DES PACKAGES ET REPERTOIRE DE TRAVAIL
################################################################################

require("raster")
### Robert J. Hijmans (2016). raster: Geographic Data Analysis and
### Modeling. R package version 2.5-8.
### https://CRAN.R-project.org/package=raster
require("rgeos")
### Roger Bivand and Colin Rundel (2017). rgeos: Interface to Geometry
### Engine - Open Source ('GEOS'). R package version 0.3-26.
### https://CRAN.R-project.org/package=rgeos
require("SDMTools")
### Jeremy VanDerWal, Lorena Falconi, Stephanie Januchowski, Luke Shoo
### and Collin Storlie (2014). SDMTools: Species Distribution Modelling
### Tools: Tools for processing data associated with species distribution
### modelling exercises. R package version 1.1-221.
### https://CRAN.R-project.org/package=SDMTools
require("rgdal")
### Roger Bivand, Tim Keitt and Barry Rowlingson (2017). rgdal: Bindings
### for the 'Geospatial' Data Abstraction Library. R package version
### 1.2-15. https://CRAN.R-project.org/package=rgdal

wd <- "/home/bioinfo/Documents/2018_M2_CamilaBenavides/SIG_QGIS/"
setwd(wd)
# list.files()

################################################################################
### [0] POINTS PARCELLES
################################################################################

### points gps des parcelles
radiusR <- seq(from = 25, to = 400, by = 25)
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
CRS_new <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
d_ch1903 <- spTransform(d, CRS_new)

################################################################################
### [1] RASTER DEPUIS QGIS 
################################################################################

### rasters 1000*1000 avec les trois categories
patch05 <- raster("p05_3k.tif")
patch07 <- raster("p07_3k.tif")
patch12 <- raster("p12_3k.tif")
patch18 <- raster("p18_3k.tif")
patch38 <- raster("p38_3k.tif")
patch40 <- raster("p40_3k.tif")
patch41 <- raster("p41_3k.tif")
patch42 <- raster("p42_3k.tif")

lpatch3k <- list(patch05, patch07, patch12, patch18, patch38, patch40, patch41,
	patch42)

# ### extraction des pourcentages de chaque categorie
RS_drone_patch_percentK1 <- lapply(1:length(lpatch3k), function(myRaster){
	posPoints <- d_ch1903
		myPoint <- posPoints[myRaster]
		indi <- sapply(radiusR, function(i){
			myRaster_k <- extract(lpatch3k[[myRaster]],
			myPoint,
			buffer = i,
			fun = function(k){
				k1 <- length(k[k == 1]) / length(k)
				return(k1)
			},
			df = FALSE)
		})
	names(indi) <- paste0("perk1_", as.character(radiusR))
	return(indi)
})
RS_drone_patch_percentK1_DF <- data.frame(
	do.call(rbind, RS_drone_patch_percentK1))
rownames(RS_drone_patch_percentK1_DF) <- as.character(pointsDF$parcelle)
RS_drone_patch_percentK2 <- lapply(1:length(lpatch3k), function(myRaster){
	posPoints <- d_ch1903
		myPoint <- posPoints[myRaster]
		indi <- sapply(radiusR, function(i){
			myRaster_k <- extract(lpatch3k[[myRaster]],
			myPoint,
			buffer = i,
			fun = function(k){
				k2 <- length(k[k == 2]) / length(k)
				return(k2)
			},
			df = FALSE)
		})
	names(indi) <- paste0("perk2_", as.character(radiusR))
	return(indi)
})
RS_drone_patch_percentK2_DF <- data.frame(
	do.call(rbind, RS_drone_patch_percentK2))
rownames(RS_drone_patch_percentK2_DF) <- as.character(pointsDF$parcelle)
RS_drone_patch_percentK3 <- lapply(1:length(lpatch3k), function(myRaster){
	posPoints <- d_ch1903
		myPoint <- posPoints[myRaster]
		indi <- sapply(radiusR, function(i){
			myRaster_k <- extract(lpatch3k[[myRaster]],
			myPoint,
			buffer = i,
			fun = function(k){
				k3 <- length(k[k == 3]) / length(k)
				return(k3)
			},
			df = FALSE)
		})
	names(indi) <- paste0("perk3_", as.character(radiusR))
	return(indi)
})
RS_drone_patch_percentK3_DF <- data.frame(
	do.call(rbind, RS_drone_patch_percentK3))
rownames(RS_drone_patch_percentK3_DF) <- as.character(pointsDF$parcelle)

################################################################################
### [2] SHAPEFILE DE CHAQUE CATHEGORIE
################################################################################

### adapted from https://rossijeanpierre.wordpress.com 
### and https://conservationecology.wordpress.com/

myShps <- list(
	c(
		"P05_k1_patchBasedClassif.shp",
		"P05_k2_patchBasedClassif.shp",
		"P05_k3_patchBasedClassif.shp"
	), 
	c(
		"P07_k1_patchBasedClassif.shp",
		"P07_k2_patchBasedClassif.shp",
		"P07_k3_patchBasedClassif.shp"
	),
	c(
		"P12_k1_patchBasedClassif.shp",
		"P12_k2_patchBasedClassif.shp",
		"P12_k3_patchBasedClassif.shp"
	),
	c(
		"P18_k1_patchBasedClassif.shp",
		"P18_k2_patchBasedClassif.shp",
		"P18_k3_patchBasedClassif.shp"
	),
	c(
		"P38_k1_patchBasedClassif.shp",
		"P38_k2_patchBasedClassif.shp",
		"P38_k3_patchBasedClassif.shp"
	),
	c(
		"P40_k1_patchBasedClassif.shp",
		"P40_k2_patchBasedClassif.shp",
		"P40_k3_patchBasedClassif.shp"
	),
	c(
		"P41_k1_patchBasedClassif.shp",
		"P41_k2_patchBasedClassif.shp",
		"P41_k3_patchBasedClassif.shp"
	),
	c(
		"P42_k1_patchBasedClassif.shp",
		"P42_k2_patchBasedClassif.shp",
		"P42_k3_patchBasedClassif.shp"
	)
)
myCat <- c("k1", "k2", "k3")
myRadiusT <- radiusR * 2 # diametre
lx <- sapply(1:length(myShps), function(numParcel){
	px <- sapply(1:length(myCat), function(numK){
		kx <- sapply(myRadiusT, function(myRadiusTi){
			teow <- readOGR(myShps[[numParcel]][numK])
			teow_utm <- spTransform(teow, CRS("+init=epsg:21037"))
			teow_utm_rst <- raster(extent(teow_utm), crs = projection(teow_utm), 
				ncol = 100, nrow = 100) #resolution
			rp <- rasterize(x = teow_utm, y = teow_utm_rst, field = 'idALL', 
				update = TRUE) #rp: le nouveau raster es cree dans le raster "VIDE": teow_utm_rst
			cpoint <- c(mean(extent(teow_utm_rst)[1:2]), 
				mean(extent(teow_utm_rst)[3:4]))
			sp <- SpatialPoints(rbind(cpoint))
			buf <- gBuffer(sp, width = myRadiusTi, quadsegs = 50)
			raster2 <- crop(rp, buf, snap = "out") 
			crop <- setValues(raster2, NA) 
			bufr <- rasterize(buf, crop) 
			out <- mask(x = raster2, mask = bufr) 
			# plot(out)
			fragInd <- ClassStat(mat = out)
			fragInd <- fragInd[fragInd$class != 0,]

			fragIndV <- list()
			fragIndV$n.patches <- sum(fragInd[,2], na.rm = TRUE)
			fragIndV$total.area <- sum(fragInd[,3], na.rm = TRUE)
			fragIndV$prop.landscape <- mean(fragInd[,4], na.rm = TRUE)
			fragIndV$patch.density <- sum(fragInd[,5], na.rm = TRUE)
			fragIndV$total.edge <- mean(fragInd[,6], na.rm = TRUE)
			fragIndV$edge.density <- mean(fragInd[,7], na.rm = TRUE)
			fragIndV$landscape.shape.index <- mean(fragInd[,8], na.rm = TRUE)
			fragIndV$largest.patch.index <- max(fragInd[,9], na.rm = TRUE)
			fragIndV$patch.area.mean <- mean(fragInd[,10], na.rm = TRUE)
			fragIndV$patch.area.sd <- sd(fragInd[,10], na.rm = TRUE)
			fragIndV$patch.area.min <- min(fragInd[,10], na.rm = TRUE)
			fragIndV$patch.area.max <- max(fragInd[,10], na.rm = TRUE)
			fragIndV$perimeter.area.frac.dim <- mean(fragInd[,14], na.rm = TRUE)
			fragIndV$perim.area.ratio.mean <- mean(fragInd[,15], na.rm = TRUE)
			fragIndV$perim.area.ratio.sd <- sd(fragInd[,15], na.rm = TRUE)
			fragIndV$perim.area.ratio.min <- min(fragInd[,15], na.rm = TRUE)
			fragIndV$perim.area.ratio.max <- max(fragInd[,15], na.rm = TRUE)
			fragIndV$shape.index.mean <- mean(fragInd[,19], na.rm = TRUE)
			fragIndV$shape.index.sd <- sd(fragInd[,19], na.rm = TRUE)
			fragIndV$shape.index.min <- min(fragInd[,19], na.rm = TRUE)
			fragIndV$shape.index.max <- max(fragInd[,19], na.rm = TRUE)
			fragIndV$frac.dim.index.mean <- mean(fragInd[,23], na.rm = TRUE)
			fragIndV$frac.dim.index.sd <- sd(fragInd[,23], na.rm = TRUE)
			fragIndV$frac.dim.index.min <- min(fragInd[,23], na.rm = TRUE)
			fragIndV$frac.dim.index.max <- max(fragInd[,23], na.rm = TRUE)
			fragIndV$total.core.area <- sum(fragInd[,27], na.rm = TRUE)
			fragIndV$prop.landscape.core <- mean(fragInd[,28], na.rm = TRUE)
			fragIndV$patch.core.area.mean <- mean(fragInd[,29], na.rm = TRUE)
			fragIndV$patch.core.area.sd <- sd(fragInd[,29], na.rm = TRUE)
			fragIndV$patch.core.area.min <- min(fragInd[,29], na.rm = TRUE)
			fragIndV$patch.core.area.max <- max(fragInd[,29], na.rm = TRUE)
			fragIndV$prop.like.adjacencies <- mean(fragInd[,33], na.rm = TRUE)
			fragIndV$aggregation.index <- mean(fragInd[,34], na.rm = TRUE)
			fragIndV$lanscape.division.index <- mean(fragInd[,35], na.rm = TRUE)
			fragIndV$splitting.index <- mean(fragInd[,36], na.rm = TRUE)
			fragIndV$effective.mesh.size <- mean(fragInd[,37], na.rm = TRUE)
			fragIndV$patch.cohesion.index <- mean(fragInd[,38], na.rm = TRUE)
			return(unlist(fragIndV))
		})
		kxx <- as.vector(kx)
		radiusNames <- myRadiusT/2
		radiusNames[radiusNames < 100] <- paste0("0", 
			radiusNames[radiusNames < 100])
		names(kxx) <- paste0(rep(rownames(kx), length(myRadiusT)), "_", 
			rep(radiusNames, each = 37))
		return(kxx)
	})
	pxx <- as.vector(px)
	names(pxx) <- paste0(rep(myCat, each = 37 * length(myRadiusT)), "_", 
		rep(rownames(px), length(myCat)))
	return(pxx)
})
print(lx)

save(lx, file = "lx.RData")

load("lx.RData")

colnames(lx) <- as.character(pointsDF$parcelle)
patchMetrics <- t(lx)
## indices par ordre alaphabetique par rayon :
## k1_n.patches_025 ; k1_n.patches_050 ; ... ; k1_n.patches_400
patchMetricsO <- patchMetrics[, order(colnames(patchMetrics))]
patchMetrics <- patchMetricsO 
## ajout de la colonne parcela avec nom des parcelles
patchMetrics <- data.frame(parcela = rownames(patchMetrics), patchMetrics)
## fusion des indices avec les variables a expliquer
phytophagy <- read.table("Tabla_integrada.txt", header = TRUE, sep = "\t")
phytophagy <- phytophagy[,c("parcela", "PhojasD", "PhojasE", "ShojasD", "IDA")]
patchMetDF <- merge(phytophagy, patchMetrics, by.x = 1, by.y = 1, all = TRUE)
## preparation des donnees pour le calcul de la matrice de correlation
varDF <- patchMetDF
explVars <- list(varDF$PhojasD, varDF$PhojasE, varDF$ShojasD, varDF$IDA)
varNames <- c("PhojasD", "PhojasE", "ShojasD", "IDA")
ncolB <- 6
ncolE <- 1781
radio <- seq(25, 400, 25)
nVar <- (ncolE - ncolB + 1)/length(radio)

## nom de chaque variable par rayon par cathegorie
XXX111XXX <- paste0(rep(c("k1_", "k2_", "k3_"), each = nVar/3), 
	unique(sapply(strsplit(colnames(varDF[,ncolB:ncolE]), split = "_"), 
	"[[", 2)))

## PDF de la matrice de correlation
pdf("corrMatPerRadiusPatchMet_varExpl_100.pdf", width = 20, height = )
	for(i in 1:length(explVars)){
		layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 2), ncol = 1))
		par(mar = c(12, 4, 4, 0))
		corVarDF <- cor(x = explVars[[i]], y = varDF[, ncolB:ncolE])#, 
			#use = "pairwise.complete.obs")
		corVarDFM <- matrix(as.vector(corVarDF), ncol = nVar)
		rownames(corVarDFM) <- paste0("r", radio)
		colnames(corVarDFM) <- XXX111XXX
		rgb.palette <- colorRampPalette(
			c("blue", "cyan", "white", "orange", "red"), space = "rgb")
		yy <- sapply(apply(t(corVarDFM), MARGIN = 1, FUN = function(k){
			which(abs(k) == max(abs(k), na.rm = TRUE))
		}), "[[", 1) - 1
		xx <- 0:(nVar - 1)
		image(t(corVarDFM), col = rgb.palette(500), axes = FALSE, 
			main = varNames[i], zlim = c(-1, 1))
		points(xx/max(xx), yy/15, pch = 16)
		axis(1, labels = XXX111XXX, at = (0:(nVar - 1))/(nVar - 1), las = 2)
		axis(2, labels = radio, at = 0:(length(radio) - 1)/(length(radio) - 1), 
			las = 2)
		par(mar = c(4, 4, 0, 90))
		image(matrix(seq(from = -1, to = 1, by = 0.1), ncol = 1), 
			col = rgb.palette(500), axes = FALSE)
		axis(1, labels = seq(from = -1, to = 1, by = 0.1), 
			at = seq(from = 0, to = 1, by = 0.05))
	}
dev.off()
