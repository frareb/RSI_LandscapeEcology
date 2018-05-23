#!/usr/bin/Rscript
###############################################################################
### [III] Statistiques descriptives donnees de niveau de phyto
### et relation niveau phyto~pratiques
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

require("MASS")
### Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with
### S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0

wd <- "/home/bioinfo/Documents/2018_M2_CamilaBenavides/SIG_QGIS/"
setwd(wd)
# list.files()

###############################################################################
### [1] CHARGEMENT DES DONNEES
###############################################################################

### [1.1] Charger les fichiers avec les donnees observees
obsAvril <- read.table("FR_abril_2017.csv", sep = ";", 
	header = TRUE, colClasses = c("factor", "factor", "factor", 
	"character", "numeric", "numeric", "numeric", "factor", "factor", 
	"numeric", "numeric", "numeric", "numeric", "factor", "character"), 
	dec = ".")
obsMai <- read.table("FR_mayo_2017.csv", sep = ";", 
	header = TRUE, colClasses = c("factor", "factor", "factor", 
	"character", "numeric", "numeric", "numeric", "factor", "factor", 
	"numeric", "numeric", "numeric", "numeric", "factor", "character"), 
	dec = ".")
obsJuin <- read.table("FR_junio_2017.csv", sep = ";", 
	header = TRUE, colClasses = c("factor", "factor", "factor", 
	"character", "numeric", "numeric", "numeric", "factor", "factor", 
	"numeric", "numeric", "numeric", "numeric", "factor", "character"), 
	dec = ".")

### [1.2] Limiter les donnees a celles effectuees dans la parcelle
obsAvril <- obsAvril[grep(obsAvril$trat, pattern ="L"),]
obsMai <- obsMai[grep(obsMai$trat, pattern ="L"),]
obsJuin <- obsJuin[grep(obsJuin$trat, pattern ="L"),]

### [1.3] limiter les donnees aux parcelles avec photo drone
parcelles <- c("P5", "P7", "P12", "P18", "P38", "P40", "P41", "P42")
obsAvril <- obsAvril[as.character(obsAvril$parcela) %in% parcelles,]
obsMai <- obsMai[as.character(obsMai$parcela) %in% parcelles,]
obsJuin <- obsJuin[as.character(obsJuin$parcela) %in% parcelles,]

### [1.4] grouper les donnees dans un objet
obsALL <- list(obsAvril, obsMai, obsJuin)
## description de l objet
sapply(obsALL, function(i){length(i[,1])}) # nombre d observations / date
sapply(obsALL, function(i){table(as.character(i$parcela))}) # nbr obs / date / parcelle
table(as.character(rbind(obsALL[[1]], obsALL[[2]], obsALL[[3]])$parcela)) # /parc

### [1.5] variables calculees ajoutees
obsALLV <- lapply(seq_along(obsALL), function(i){
	obsALL[[i]] <- cbind(obsALL[[i]], 
		PhojasD = obsALL[[i]]$hojas_dañadas / obsALL[[i]]$hojas_total * 100)
	obsALL[[i]] <- cbind(obsALL[[i]], 
		PhojasE = obsALL[[i]]$hojas_enrolladas / obsALL[[i]]$hojas_tot_enr * 100)
	obsALL[[i]] <- cbind(obsALL[[i]], 
		ShojasD = (obsALL[[i]]$hojas_dañadas / obsALL[[i]]$hojas_total) * 
		as.numeric(obsALL[[i]]$superficie_daño))
	obsALL[[i]] <- cbind(obsALL[[i]], 
		IDA = (obsALL[[i]]$num_ramificaciones / obsALL[[i]]$tamaño_tallo))
})

obsALLVm_itk2 <- read.table("obsALLVm_itk.txt", header = TRUE, sep = "\t")

###############################################################################
### [2] STATISTIQUES DESCRIPTIVES
###############################################################################

#### 2.1. boxplot indices dommage PER mois d'echantillonage
###############################################################################
par(mfrow = c(1, 4))
for(j in c(16:19)){
  boxplot(lapply(obsALLV, function(i){
    i[,j]
  }), 
  names = c("Avril", "Mai", "Juin"), ylab = names(obsALLV[[1]])[j]
  )
}

#### 2.2. boxplot indices dommage PER parcelle
###############################################################################
obsALLV[[1]] <- cbind(obsALLV[[1]], mes = "Avril")
obsALLV[[2]] <- cbind(obsALLV[[2]], mes = "Mai")
obsALLV[[3]] <- cbind(obsALLV[[3]], mes = "Juin")
obsALLVm <- rbind(obsALLV[[1]], obsALLV[[2]], obsALLV[[3]])
dfCol <- data.frame(myCol = rainbow(8), parcelles)

par(mfrow = c(1, 4))
for(j in c(16:19)){
	splitObsALLVm <- split(obsALLVm[,j], as.character(obsALLVm$parcela))
	splitObsALLVm <- splitObsALLVm[order(sapply(splitObsALLVm, mean, 
		na.rm = TRUE))]
	boxplot(splitObsALLVm, col = as.character(merge(names(splitObsALLVm), 
		dfCol, by.x = 1, by.y = 2, sort = FALSE)[,2]), 
		ylab = names(obsALLVm)[j])
}

### 2.3. Statistiques descriptives de phyto 
###############################################################################
for(i in 16:19){
	print(names(obsALLVm)[i])
	print(c(summary(obsALLVm[,i]), sd = sd(obsALLVm[,i], na.rm = TRUE))) # OJO 
	# VALEUR 6 DE SHD!!!!
}

###############################################################################
### [3] CHARGEMENT DONNES ITK
###############################################################################

### [3.1] Donnees colectees
###############################################################################
itk <- data.frame(
	parcelle = c("P5", "P7", "P12", "P18", "P38", "P40", "P41", "P42"), 
	longitude = c(-78.78509333, -78.80273167, -78.79859667, -78.63821333, 
		-78.588908, -78.743102, -78.730944, -78.76065), 
	latitude = c(-2.108118333, -2.088856667, -2.100291667, -0.963956667, 
		-0.907498, -0.824985, -0.815995, -2.008947), 
	prec1 = c("Cebada", "Descanso", "Chocho", "Descanso", "Maiz", "Descanso", 
		NA, "Cebada"),
	prec2 = c("Chocho", "Descanso", "Descanso", "Descanso", "Chocho", 
		"Descanso", NA, "Quinua"),
	phyto = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, NA, TRUE), 
	varie = c("Andino", "Andino", "Andino", "Andino", "Andino", "Andino", 
		"Andino", "Andino"), 
	dateSem = c("17/02/2017", "23/02/2017", "21/02/2017", "09/03/2017", 
		"23/02/2017", "19/02/2017", "22/02/2017", NA),
	provincia = c("Chimborazo", "Chimborazo", "Chimborazo", "Cotopaxi", 
		"Cotopaxi", "Cotopaxi", "Cotopaxi", "Chimborazo"),
	dateSem_mois = c("Fevrier", "Fevrier", "Fevrier", "Mars", "Fevrier", 
		"Fevrier", "Fevrier", NA)
)

### [3.2] analyses preliminaires
###############################################################################
obsALLVm_itk <- merge(obsALLVm, itk, by.x = 1, by.y = 1)
obsALLVm_itk$dateSem <- as.Date(obsALLVm_itk$dateSem, format = "%d/%m/%Y")
obsALLVm_itk$dateSemDif <- obsALLVm_itk$dateSem - min(obsALLVm_itk$dateSem, na.rm = TRUE)
