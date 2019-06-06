rm(list=ls())
setwd("~/R/numerico2")
library(pixmap)
library(MASS)
files <- list.files(path=".", pattern=".pgm",all.files=T, full.names=F, no.. = T) 
list_of_images = lapply(files, read.pnm)

matriz2<-cbind(c(list_of_images[[1]]@grey),c(list_of_images[[2]]@grey),c(list_of_images[[3]]@grey),c(list_of_images[[4]]@grey),c(list_of_images[[5]]@grey),
               c(list_of_images[[6]]@grey),c(list_of_images[[7]]@grey),c(list_of_images[[8]]@grey),c(list_of_images[[9]]@grey),c(list_of_images[[10]]@grey),
               c(list_of_images[[11]]@grey))

colnames(matriz2)<-c("sujeto1","sujeto2","sujeto3","sujeto4","sujeto5","sujeto6","sujeto7","sujeto8",
                     "sujeto9","sujeto10","sujeto11")

matriz3<-cbind(c(list_of_images[[12]]@grey),c(list_of_images[[13]]@grey),c(list_of_images[[14]]@grey),c(list_of_images[[15]]@grey),c(list_of_images[[16]]@grey),
               c(list_of_images[[17]]@grey),c(list_of_images[[18]]@grey),c(list_of_images[[19]]@grey),c(list_of_images[[20]]@grey),c(list_of_images[[21]]@grey),
               c(list_of_images[[22]]@grey),c(list_of_images[[23]]@grey),c(list_of_images[[24]]@grey))

etiquetas_test<-c(1,4,5,8,10,4,5,6,7,8,3,10,11)

setwd("~/ciencia de datos/Proyecto final")
library(Rcpp)
library(RcppArmadillo)
library(devtools)
sourceCpp("costos.cpp")
source("distancias.R")

LNMF2<-function(V,rango,tol=0.001, max.iter=10,semilla=100){
  set.seed(semilla)
  W<-matrix(sample(seq(1,255,1),dim(V)[1]*rango,replace=TRUE), ncol=rango)
  H<-matrix(sample(seq(1,255,1),dim(V)[2]*rango,replace=TRUE), nrow=rango)
  costoInicial<-costoNMF(V,W,H)
  #print(costoInicial)
  for (i in 1:max.iter){
    resultado=updateNMF(V,W,H)
    #print(costoLNMF(V, resultado[[1]],resultado[[2]]))
    if (abs(costoNMF(V, resultado[[1]],resultado[[2]])-costoInicial)<tol){
      break;
    }
    costoInicial<-costoNMF(V, resultado[[1]],resultado[[2]])
  }
  resultado
}





predictNMF<-function(modelo, train, test, distancia="euclid"){
  require(MASS)
  n<-dim(test)[2]
  clasif<-rep(0,n)
  media<-apply(train,1,mean)
  baseMoore<-ginv(modelo[[1]])
  train.red<-baseMoore%*%(train)
  test.red<-baseMoore%*%(test)
  if (distancia=="euclid"){
    for (i in 1:n){
      clasif[i]<-distancia(test.red[,i], train.red)
    }
  }
  else{
    for (i in 1:n){
      clasif[i]<-distCor(test.red[,i], train.red)
    }
  }
  clasif
}

predictLNMF<-function(modelo, train, test, distancia="euclid"){
  require(MASS)
  n<-dim(test)[2]
  clasif<-rep(0,n)
  media<-apply(train,1,mean)
  baseMoore<-ginv(modelo[[1]])
  train.red<-baseMoore%*%(train)
  test.red<-baseMoore%*%(test)
  if (distancia=="euclid"){
    for (i in 1:n){
      clasif[i]<-distancia(test.red[,i], train.red)
    }
  }
  else{
    for (i in 1:n){
      clasif[i]<-distCor(test.red[,i], train.red)
    }
  }
  clasif
}

for (i in 2:11){
  a<-LNMF2(matriz2,i,100)
  b<-predictLNMF(a,matriz2,matriz3,distancia = "corr")
  c<-sum(ifelse(etiquetas_test==b,1,0))/13
  print(paste(i," ",c*100))
}

library(imager)
plot(as.cimg(matrix((a[[1]]%*%a[[2]])[,5],nrow=112)))
plot(as.cimg(matrix(a[[1]][,1],ncol=112, nrow=92)))
########################## Frutas
rm(list=ls())
setwd("~/R/numerico")
library(imager)
files <- list.files(path=".", pattern=".jpg",all.files=T, full.names=F, no.. = T) 
list_of_images = lapply(files, load.image)

matriz<-cbind(c(grayscale(list_of_images[[1]])),c(grayscale(list_of_images[[2]])),c(grayscale(list_of_images[[3]])),c(grayscale(list_of_images[[4]])),c(grayscale(list_of_images[[5]])),c(grayscale(list_of_images[[6]])),
              c(grayscale(list_of_images[[7]])),c(grayscale(list_of_images[[8]])),c(grayscale(list_of_images[[9]])), c(grayscale(list_of_images[[10]])),c(grayscale(list_of_images[[11]])),
              c(grayscale(list_of_images[[12]])),c(grayscale(list_of_images[[13]])),c(grayscale(list_of_images[[14]])), c(grayscale(list_of_images[[15]])),c(grayscale(list_of_images[[16]])),
              c(grayscale(list_of_images[[17]])),c(grayscale(list_of_images[[18]])),c(grayscale(list_of_images[[19]])), c(grayscale(list_of_images[[20]])),c(grayscale(list_of_images[[21]])))

colnames(matriz)<-c("Braeburn","Golden","Albaricoque","Aguacate","Carambolo","Arándano","Kiwi","Naranja",
                    "Durazno","Piña","Fresa","Durazno","Aguacate","Carambolo","Kiwi", "Piña", "Kiwi", "Naranja", 
                    "Braeburn","Arándano", "Carambolo")

etiquetas_test=c(10, 4,5,7,11,7,8,1,6,5)
for (i in 2:11){
  a<-NMF2(matriz[,1:11],i,100)
  b<-predictNMF(a,matriz[,1:11],matriz[,12:21],distancia = "corr")
  c<-sum(ifelse(etiquetas_test==b,1,0))/10
  print(paste(i," ",c))
}
