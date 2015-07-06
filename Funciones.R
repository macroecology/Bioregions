### FUNCIONES ##################################

### Generar matriz nodos filogenéticos x taxa 
get_taxa_nodes<- function (matrix, tree){  
  nt<-length(tree$tip.label) # número de taxa 28
  nn <- tree$Nnode # número de internal nodes 26
  ncol(matrix) # número de taxa en matriz 54
  nni<-c(1:(nt+nn)) # id nodes 1:54
  names <- tree$tip.label # nombre de los taxones
  bucle <- function(x){
    a <- tips(tree, x) # tips() da de una filog los descendientes de un nodo determinado
    b <- as.numeric(names %in% a) 
  }
  res <- sapply(as.list(nni), bucle, simplify = "matrix")
      # recorre todos los nodos de la filogenia 
      # primero ascendentemente en los tips (nodos terminales, especies) y luego
      # jerárquicamente en la fil. (de más antiguo a más reciente)
      # también de abajo a arriba,                              
      # dando un valor 1 a los descendientes de cada uno de esos nodos,
      # pero sólo para los nodos internos.
  colnames(res) <- nni # nombre de nodos en columnas
  rownames(res) <- names # nombre de taxa en filas
  res # matriz nodos x taxa
}


### Generar matriz presencia/ausencia de nodos flogenéticos, celdas x nodos
get_nodes_cells <- function(res, matrix, names){ # res es matriz nodos x taxa, matrix 
  # la matriz de presencia/ausencia (celdas geográficas x taxa) y names el nombre de los
  # taxones en la filogenia
names2 <- colnames(matrix) # nombre de los taxones en la matriz
matrix2 <- ifelse(matrix == 1, TRUE, FALSE) # se convierte la matriz de 1/0 en matriz lógica
  
posbu <- apply(matrix2, 1, function(x){names %in% names2[x]})

rfin <- t(apply(posbu, 2, function(x){colSums(res[x, , drop = FALSE])}))
colnames(rfin) <- nni # los nodos
rownames(rfin) <- rownames(matrix) # las celdas geográficas. Así se le da más peso a los nodos antiguos.
rfin2 <- ifelse(rfin > 0, 1, 0) # y así se le da igual peso a todos los nodos. (*)

nam<-data.frame(c(rownames(rfin),paste("x",colnames(rfin),sep=""))) # nombres para meter en pajek y asi
  # generar el input para infomap. Importante añadir la x antes del número de cada nodo!
res<- list (rfin=rfin, rfin2=rfin2, nam=nam) 
}

# (*) Este sería el primer punto de partida después de entender bien cómo se generarn posbu y rfin. Habría
# que incorporar un argumento más a la función con el que se pueda elegir el tipo de pesos que queremos dar, a elegir
# entre: igual a todos los nodos, más a los más antiguos  o en función de GND. Para el último solamente sería multiplicar
# 1 / el vector "pesos" por la matriz rfin. Una vez está rfin con el peso que queremos los siguientes pasos son:

# Generar input para infomap en R, información:
 # Generar pajek formats http://igraph.org/r/
 # En http://www.mapequation.org/code.html , mirar la sección dedicada a R en el apartado 2-Runing-with other languages  

# Correr infomap en R si es posible, y si no, llamar a C++ desde R para correrlo. Información:
 # la función infomap.community del paquete igraph parece que computa infomap directamente.

# Guardar resultados en forma de tabla (asignación de módulos por nodo y también por celda geográfica) y de mapa, 
# visualizar con raster las biorregiones.

# Hacer un loop para que el proceso se repita 100 veces con la misma matriz presencia/ausencia de taxa por celdas geográficas
# pero 100 filogenias diferentes.

# Guardar resultados de todas las iteraciones en la misma tabla y en el mismo mapa (tengo que pensar cómo hacer 
# esta representación gráfica).





