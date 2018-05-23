#!/usr/bin/Rscript
###############################################################################
### [VIII] Script des analyses: 
### refuges vs. dommage
### var LANDSAT, pixel, DSM(surface) vs. dommage
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

require("raster")
### Robert J. Hijmans (2016). raster: Geographic Data Analysis and
### Modeling. R package version 2.5-8.
### https://CRAN.R-project.org/package=raster

wd <- "/home/bioinfo/Documents/2018_M2_CamilaBenavides/SIG_QGIS/"
setwd(wd)
# list.files()

###########################################################################
### [1] CORRELATION ENTRE REFUGES ET DOMMAGE
###########################################################################

refuges <- read.table("DfRefuge.txt", header = TRUE, sep = "\t")
explVars <- list(refuges$PhojasD, refuges$PhojasE, refuges$ShojasD, refuges$IDA)
varNames <- c("PhD", "PhE", "ShD", "IDA")
ncolB <- 6
ncolE <- 293
radio <- seq(25, 400, 25)
nVar <- (ncolE - ncolB + 1)/length(radio)

pdf("corrMatRefuges_varExpl.pdf", width = 5, height = 5)
	par(mar = c(4, 4, 4, 0))
	for(i in 1:length(explVars)){
		corRefuges <- cor(x = explVars[[i]], y =  refuges[, 6:ncol(refuges)])
		prof <- seq(from = 0.3, to = 2, by = 0.1)
		reso <- seq(form = 50, to = 200, by = 10)
		corRefugesM <- matrix(as.vector(corRefuges), ncol = 18, byrow = TRUE)
		rownames(corRefugesM) <- paste0("r", reso)
		colnames(corRefugesM) <- paste0("p", prof)
		corRefugesPval <- sapply(6:ncol(refuges), function(j){
			cor.test(x = explVars[[i]], y = refuges[,j])$p.value
		})
		corRefugesPvalM <- matrix(as.vector(corRefugesPval), ncol = 18, byrow = TRUE)
		rownames(corRefugesPvalM) <- paste0("r", reso)
		colnames(corRefugesPvalM) <- paste0("p", prof)
		rgb.palette <- colorRampPalette(
			c("blue", "cyan", "white", "orange", "red"), space = "rgb")
		yy <- (apply(t(corRefugesM), MARGIN = 1, FUN = function(k){
			which(abs(k) == max(abs(k)))
		}) - 1)
		xx <- 0:(nVar - 1)
		image(corRefugesM, col = rgb.palette(500), axes = FALSE, main = varNames[i], 
			zlim = c(-1, 1))
		# points(y = xx/max(xx), x = yy/15, pch = 16)
		contour(corRefugesM, add = TRUE, lwd = 2)
		axis(1, labels = paste0("r", reso), at = (0:15)/(15), las = 2)
		axis(2, labels = paste0("p", prof), at = 0:(17)/(17), las = 2)
	}
dev.off()

###########################################################################
### [2] CORRELATION ENTRE INDICES LANDSAT, PIXEL, DSM ET DOMMAGE
###########################################################################

varDF <- read.table("Tabla_integrada2.txt", header = TRUE, sep = "\t")
explVars <- list(varDF$PhojasD, varDF$PhojasE, varDF$ShojasD, varDF$IDA)
varNames <- c("PhojasD", "PhojasE", "ShojasD", "IDA")
ncolB <- 15
ncolE <- 222
radio <- seq(from = 25, to = 400, by = 25)
nVar <- (ncolE - ncolB + 1)/length(radio)
nombreCol <- c("TC_wet", "TC_green", "TC_bright",
               "NDVI", "pixNS_SD", "pixNS_MAX", "pixS_k1", 
               "pixS_k2", "pixS_k3", "DSM_hill", "DSM_rough", 
               "DSM_aspect", "DSM_slope")

pdf("corrMatPerRadiusLandDSM_varExpl2.pdf", width = 4, height = 5)
	par(mar = c(8, 4, 4, 1))
	for(i in 1:length(explVars)){
		corVarDF <- cor(x = explVars[[i]], y =  varDF[, ncolB:ncolE])
		corVarDFM <- matrix(as.vector(corVarDF), ncol = nVar)
		rownames(corVarDFM) <- paste0("r", radio)
		colnames(corVarDFM) <- nombreCol
		corVarPval <- sapply(ncolB:ncolE, function(j){
			cor.test(x = explVars[[i]], y = varDF[,j])$p.value
		})
		corVarPvalM <- matrix(as.vector(corVarPval), ncol = nVar, byrow = TRUE)
		rownames(corVarPvalM) <- paste0("r", radio)
		colnames(corVarPvalM) <- nombreCol
		rgb.palette <- colorRampPalette(
			c("blue", "cyan", "white", "orange", "red"), space = "rgb")
		yy <- (apply(t(corVarDFM), MARGIN = 1, FUN = function(k){
			which(abs(k) == max(abs(k)))
		}) - 1)
		xx <- 0:(nVar - 1)
		image(t(corVarDFM), col = rgb.palette(500), axes = FALSE, main = varNames[i], 
			zlim = c(-1, 1))
		points(xx/max(xx), yy/15, pch = 16)
		axis(1, labels = nombreCol, at = (0:(nVar - 1))/(nVar - 1), las = 2)
		axis(2, labels = radio, at = 0:(length(radio) - 1)/(length(radio) - 1), las = 2)
	}
dev.off()
