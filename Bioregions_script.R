## PRUEBAS BIORREGIONES Carnivora en Europa ##############
# Se cargan los paquetes necesarios
library(ape)
library(nlme)
library(geiger)
library(maptools)
library(sp)
library(raster)
library(maps)
library(letsR)
library(nodiv)

# Se abren los datos y se crea la matriz presencia/ausencia

tree <- read.nexus("C:/Users/Anbrial/Desktop/TalleR/TalleR_Hyenidae/Data/phylcarnivora.txt")
plot(tree)

dc <- readShapePoly("C:/Users/Anbrial/Desktop/TalleR/TalleR_Hyenidae/Data/carnivora/EU_extant_carnivora.shp") 
pc <- lets.presab(dc, xmn = -30, xmx = 75, ymn = 30, ymx = 90, resol=1,
                  crs=CRS("+proj=longlat +datum=WGS84"))
x11()
plot(pc)

pcm <- pc[[1]]
coord <- pcm[,c(1,2)]
dim(coord)

# Matching entre filogenia y matriz
colnames(pcm) <- gsub(" ", "_", colnames(pcm))
km<-data.frame(matrix(nrow=ncol(pcm[,-c(1:2)]),ncol=1)) # matriz de taxa de pm x 1 col
kf<-data.frame(matrix(ncol=1,nrow=length(tree$tip.label))) # matriz de taxa de filogenia x 1 col
names(kf)<-names(km)<-c("names") # damos un nombre unificado a los nombres de filas y columnas de kf ykm.
km[,2]<-1
kf[,2]<-1
km[,1]<-colnames(pcm[,-c(1:2)]) # llenamos la única columna de km con los nombres de las familias de la matriz pm.
kf[,1]<-tree$tip.label # llenamos la única columna de kf con los nombres de las familias de la filogenia
kd<-merge(kf,km,by="names", all.x=T,all.y=T) # unimos ambas columnas en un mismo data frame
km<-kd[is.na(kd[,3]),1] # especies presentes en la filogenia pero no en la matriz
kf<-kd[is.na(kd[,2]),1] # especies presentes en la matriz pero no en la filogenia
## Quito taxa matriz
pcm.def <- pcm[,-c(1:2)][,!colnames(pcm[,-c(1:2)])%in% kf] # sale bien pero sin las coord
## Quito taxa filogenia
tree.def <- drop.tip(tree,km)
## Miro si ya tienen el mismo número de especies:
ncol(pcm.def) # 4701 especies
length(tree.def$tip.label) # 4701 especies.

pcm.def <- cbind(coord, pcm.def)

# Análisis de distribuciones en nodos: paquete nodiv
# cargo los paquetes necesarios 
library(picante)
# creo el nodiv object
nod.ma <- nodiv_data(tree.def, pcm.def[,-c(1,2)], coord,
                     proj4string_in = CRS("+proj=longlat +datum=WGS84"), 
                     type = c("grid"), shape = NULL)

# hago el análisis
nodan.ma <- Node_analysis(nod.ma, repeats = 100, method = c("rdtable", "quasiswap"), cores = 1, log_parallel_progress = FALSE)
str(nodan.ma)
x11()
plot(nodan.ma)

pesos<- GND(nodan.ma, node = NULL) 


#### MAP EQUATION: PREPARANDO LOS DATOS ################################
# Se cargan los paquetes necesarios
library(ape)
library(geiger)

res0 <- get_taxa_nodes (pcm.def, tree.def)

res1 <- get_nodes_cells(res0, pcm.def, tree$tip.label)
# funciones detalladas en el script Funciones.R

str(res1)
head(res1$rfin2)
dim(res1$rfin)

# Siguientes pasos detallados en el script "Funciones"


