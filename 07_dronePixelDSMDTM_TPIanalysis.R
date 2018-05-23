#!/usr/bin/Rscript
###########################################################################
### [VII] DSM DTM : calcul d un nouvel indice de REFUGES
### => change resolution des images puis exportation des donnees 
### sur des rayons de 25 a 400 m puis calcul de l indice
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
### [1] calculer le TPI des DSM par resolution et profondeur
###########################################################################
myRasterP5_DSM <- raster("DRONE/P5_dsm.tif")
myRasterP7_DSM <- raster("DRONE/P7_dsm.tif")
myRasterP12_DSM <- raster("DRONE/P12_dsm.tif")
myRasterP18_DSM <- raster("DRONE/P18_dsm.tif")
myRasterP38_DSM <- raster("DRONE/P38_dsm.tif")
myRasterP40_DSM <- raster("DRONE/P40_dsm.tif")
myRasterP41_DSM <- raster("DRONE/P41_dsm.tif")
myRasterP42_DSM <- raster("DRONE/P42_dsm.tif")

listRaster <- list(myRasterP5_DSM, myRasterP7_DSM, myRasterP12_DSM, 
                   myRasterP18_DSM, myRasterP38_DSM, myRasterP40_DSM, 
				   myRasterP41_DSM, myRasterP42_DSM)

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

prof <- seq(from = 0.2, to = 4, by = 0.1) # profondeur
reso <- seq(from = 10, to = 500, by = 10) # rang

DfRefuge <- sapply(seq_along(listRaster), function(i){
  xxx <- sapply(reso, function(myReso){
    myRasterPX_DSM_fX <- aggregate(listRaster[[i]], fact = myReso) #1 
    TPIX <- terrain(x = myRasterPX_DSM_fX, opt = "TPI")
    myPoint <- d_ch1903[i] #1
    myProfs <- sapply(prof, function(myProf){
      TPIXb <- TPIX
      TPIXb[TPIXb >  - myProf] <- 0
      TPIXb[TPIXb <= - myProf] <- 1
      # plot(TPIXb)
      return(mean(unlist(extract(TPIXb, myPoint, buffer = 300, df = FALSE))))
    })
    return(myProfs)
  })
  profReso <- as.vector(xxx)
  names(profReso) <- paste0("reso_",rep(reso, each = length(prof)), "_prof_", prof)
  return(profReso)
})

DfRefuge
colnames(DfRefuge) <- as.character(pointsDF$parcelle)
save(DfRefuge, file = "DfRefuge.RData")


mRefuge <- matrix(DfRefuge[,1], ncol = length(reso), byrow = FALSE)
colnames(mRefuge) <- reso
rownames(mRefuge) <- prof
image(mRefuge, col = terrain.colors(50))
contour(mRefuge, add = TRUE)
boxplot(DfRefuge)
