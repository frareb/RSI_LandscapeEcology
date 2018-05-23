#!/usr/bin/Rscript
###############################################################################
### [I] Script pour le pre traitement des images satellite LANDSAT 8
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
### [1] PRE-TRAITEMENT DES IMAGES : Indices spectraux et Tasseled Cap
###############################################################################

### [1.0] Chargement des bandes LANDSAT 8 (bandes 2 a 7)
###############################################################################
myRasBrick <- brick(x = list(
	raster("LC08_L1TP_010061_20170920_20171012_01_T1_B2.TIF"), 
	raster("LC08_L1TP_010061_20170920_20171012_01_T1_B3.TIF"),
	raster("LC08_L1TP_010061_20170920_20171012_01_T1_B4.TIF"),
	raster("LC08_L1TP_010061_20170920_20171012_01_T1_B5.TIF"),
	raster("LC08_L1TP_010061_20170920_20171012_01_T1_B6.TIF"),
	raster("LC08_L1TP_010061_20170920_20171012_01_T1_B7.TIF")
	))
	
### [1.1] Obtenir les indices "tasseled caps"
###############################################################################
myTassCap <- tasseledCap(img = myRasBrick, sat = "Landsat8OLI")
plot(myTassCap)
TC_C2 <- writeRaster(myTassCap, 
	filename = "TC_C.tif", 
	format = "GTiff", 
	bylayer = TRUE,
	suffix = 'numbers', 
	bandorder = 'BIL', 
	overwrite = TRUE)
	
### [1.2] Obtenir les indices spectraux (NDVI, ...)
###############################################################################
mySpecInd <- spectralIndices(img = myRasBrick, blue = 1, green = 2, red = 3, 
	nir = 4, swir1 = NULL, swir2 = 5, swir3 = 6)
plot(mySpecInd)
TC_C3 <- writeRaster(mySpecInd, 
	filename = "TC_C.tif", 
	format = "GTiff", 
	bylayer = TRUE,
	suffix = 'numbers', 
	bandorder = 'BIL', 
	overwrite = TRUE)
