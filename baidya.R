rm(list=ls())
setwd("~/ciencia de datos/Proyecto final")
library(Rcpp)
library(RcppArmadillo)
library(devtools)
library(MASS)
sourceCpp("costos.cpp")
source("distancias.R")
source("prediccion.R")
etiquetas <- read.csv("~/ciencia de datos/Tarea4/class_img_exp.dat", sep="")
colnames(etiquetas)<-c("files","semanticExp","fileExp")

library(tiff)
library(imager)
setwd("~/ciencia de datos/Tarea4/img_expression/first/train")
#Leyendo imagenes
path.train="~/ciencia de datos/Tarea4/img_expression/first/train"
path.test="~/ciencia de datos/Tarea4/img_expression/first/test"
files_train <- list.files(path=path.train, pattern=".tiff",all.files=T, full.names=F, no.. = T) 
files_test <- list.files(path=path.test, pattern=".tiff",all.files=T, full.names=F, no.. = T) 

images_train = lapply(files_train, readTIFF)
setwd("~/ciencia de datos/Tarea4/img_expression/first/test")
images_test = lapply(files_test, readTIFF)

#Generamos las etiquetas
etiquetas_train<- as.data.frame(files_train)
colnames(etiquetas_train)<-"files"
etiquetas_train <- merge(etiquetas_train,etiquetas,by="files")

etiquetas_test<- as.data.frame(files_test)
colnames(etiquetas_test)<-"files"
etiquetas_test <- merge(etiquetas_test,etiquetas,by="files")

#Transformando las imagenes a un formato imager y en escala de grises
for (i in 1:length(images_train)){
  images_train[[i]]<-grayscale(as.cimg(images_train[[i]]))
}


for (i in 1:length(images_test)){
  images_test[[i]]<-grayscale(as.cimg(images_test[[i]]))
}
#Generando la matriz de expresiones
matriz_expres<-matrix(0L,ncol=length(images_train),nrow=dim(images_train[[1]])[1]*dim(images_train[[2]])[1])
for (i in 1:length(images_train)){
  matriz_expres[,i]<-c(images_train[[i]])
}

#Generando la matriz de expresiones de prueba
matriz_expres_test<-matrix(0L,ncol=length(images_test),nrow=dim(images_test[[1]])[1]*dim(images_test[[2]])[1])
for (i in 1:length(images_test)){
  matriz_expres_test[,i]<-c(images_test[[i]])
}

###Función que genera la descomposión NMF
nonNegative<-function(V,rango,tol=0.001, max.iter=10,semilla=100){
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

####Función que genera la descomposión local NMF
LnonNegative<-function(V,rango,tol=0.001, max.iter=10,semilla=100){
  set.seed(semilla)
  W<-matrix(sample(seq(1,255,1),dim(V)[1]*rango,replace=TRUE), ncol=rango)
  H<-matrix(sample(seq(1,255,1),dim(V)[2]*rango,replace=TRUE), nrow=rango)
  costoInicial<-costoLNMF(V,W,H)
  #print(costoInicial)
  for (i in 1:max.iter){
    resultado=updateLNMF(V,W,H)
    #print(costoLNMF(V, resultado[[1]],resultado[[2]]))
    if (abs(costoLNMF(V, resultado[[1]],resultado[[2]])-costoInicial)<tol){
      break;
    }
    costoInicial<-costoLNMF(V, resultado[[1]],resultado[[2]])
  }
  resultado
}

#Funcion que genera las predicciones tanto para NMF como para LNMF con base en alguna
#de las distancias 
predecir<-function(modelo, train, test, distancia="euclid"){
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
  else if(distancia=="corr"){
    for (i in 1:n){
      clasif[i]<-distCor(test.red[,i], train.red)
    }
  }
  else{
    print("No es una función de distancia valida") 
  }
  clasif
}

grafica<-matrix(0L, ncol=8,nrow=2)
a<-rep(0,dim(matriz_expres)[2])
b<-rep(0,dim(matriz_expres)[2])
contador<-1
for (i in c(10,15,20,25,30,35,40,45)){
  NMFexpr<-nonNegative(matriz_expres,i,100)
  LNMFexpr<-LnonNegative(matriz_expres,i,100)
  
  clasifNMF<-predecir(NMFexpr, matriz_expres, matriz_expres,distancia = "corr")
  clasifLNMF<-predecir(LNMFexpr, matriz_expres, matriz_expres, distancia = "corr")
  for (j in 1:length(clasifNMF)){
    a[j]<-etiquetas_train[clasifNMF[j],2]
    b[j]<-etiquetas_train[clasifLNMF[j],2]
  }
  grafica[1, contador]<-sum(a==as.numeric(as.factor(etiquetas_train[,2])))/length(a)
  grafica[2, contador]<-sum(b==as.numeric(as.factor(etiquetas_train[,2])))/length(b)
  contador<-contador+1
}

plot(c(10,15,20,25,30,35,40,45),grafica[1,], type="l", main="Train set", sub="Clasificación semantic", col=rgb(0.2,0.4,0.1,0.7) , lwd=3 , pch=17,xlab="Rango",ylab="Nivel de precisión",ylim=c(min(grafica[1,],grafica[2,]),max(grafica[1,],grafica[2,])))
lines(c(10,15,20,25,30,35,40,45),grafica[2,],col=rgb(0.9,0.4,0.1,0.7) , lwd=3 , pch=19 , type="b")
legend("topleft", 
       legend = c("NMF", "Local NMF"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.9,0.4,0.1,0.7)), 
       pch = c(17,19,19), 
       bty = "n", 
       pt.cex = 1, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.001, 0.001))

grafica2<-matrix(0L, ncol=8,nrow=2)
a<-rep(0,dim(matriz_expres)[2])
b<-rep(0,dim(matriz_expres)[2])
contador<-1
for (i in c(10,15,20,25,30,35,40,45)){
  NMFexpr<-nonNegative(matriz_expres,i,1000)
  LNMFexpr<-LnonNegative(matriz_expres,i,1000)
  
  clasifNMF<-predecir(NMFexpr, matriz_expres, matriz_expres_test,distancia = "corr")
  clasifLNMF<-predecir(LNMFexpr, matriz_expres, matriz_expres_test,distancia = "corr")
  for (j in 1:length(clasifNMF)){
    a[j]<-etiquetas_train[clasifNMF[j],2]
    b[j]<-etiquetas_train[clasifLNMF[j],2]
  }
  grafica2[1, contador]<-sum(a==as.numeric(as.factor(etiquetas_test[,2])))/length(a)
  grafica2[2, contador]<-sum(b==as.numeric(as.factor(etiquetas_test[,2])))/length(b)
  contador<-contador+1
}

plot(c(10,15,20,25,30,35,40,45),grafica2[1,], type="l", main="Test set", sub="Clasificación semantic", col=rgb(0.2,0.4,0.1,0.7) , lwd=3 , pch=17,xlab="Rango",ylab="Nivel de precisión",ylim=c(min(grafica2[1,],grafica2[2,]),max(grafica2[1,],grafica2[2,])))
lines(c(10,15,20,25,30,35,40,45),grafica2[2,],col=rgb(0.9,0.4,0.1,0.7) , lwd=3 , pch=19 , type="b")
legend("topleft", 
       legend = c("NMF", "Local NMF"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.9,0.4,0.1,0.7)), 
       pch = c(17,19,19), 
       bty = "n", 
       pt.cex = 1, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.001, 0.001))


#### Segunda categoria
grafica3<-matrix(0L, ncol=8,nrow=2)
a<-rep(0,dim(matriz_expres)[2])
b<-rep(0,dim(matriz_expres)[2])
contador<-1
for (i in c(10,15,20,25,30,35,40,45)){
  NMFexpr<-nonNegative(matriz_expres,i,100)
  LNMFexpr<-LnonNegative(matriz_expres,i,100)
  
  clasifNMF<-predecir(NMFexpr, matriz_expres, matriz_expres,distancia = "corr")
  clasifLNMF<-predecir(LNMFexpr, matriz_expres, matriz_expres, distancia = "corr")
  for (j in 1:length(clasifNMF)){
    a[j]<-etiquetas_train[clasifNMF[j],3]
    b[j]<-etiquetas_train[clasifLNMF[j],3]
  }
  grafica3[1, contador]<-sum(a==as.numeric(as.factor(etiquetas_train[,3])))/length(a)
  grafica3[2, contador]<-sum(b==as.numeric(as.factor(etiquetas_train[,3])))/length(b)
  contador<-contador+1
}

plot(c(10,15,20,25,30,35,40,45),grafica3[1,], type="l", main="Train set", sub="Clasificación file", col=rgb(0.2,0.4,0.1,0.7) , lwd=3 , pch=17,xlab="Rango",ylab="Nivel de precisión",ylim=c(min(grafica3[1,],grafica3[2,]),max(grafica3[1,],grafica3[2,])))
lines(c(10,15,20,25,30,35,40,45),grafica3[2,],col=rgb(0.9,0.4,0.1,0.7) , lwd=3 , pch=19 , type="b")
legend("topleft", 
       legend = c("NMF", "Local NMF"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.9,0.4,0.1,0.7)), 
       pch = c(17,19,19), 
       bty = "n", 
       pt.cex = 1, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.001, 0.001))

grafica4<-matrix(0L, ncol=8,nrow=2)
a<-rep(0,dim(matriz_expres)[2])
b<-rep(0,dim(matriz_expres)[2])
contador<-1
for (i in c(10,15,20,25,30,35,40,45)){
  NMFexpr<-nonNegative(matriz_expres,i,1000)
  LNMFexpr<-LnonNegative(matriz_expres,i,1000)
  
  clasifNMF<-predecir(NMFexpr, matriz_expres, matriz_expres_test,distancia = "corr")
  clasifLNMF<-predecir(LNMFexpr, matriz_expres, matriz_expres_test,distancia = "corr")
  for (j in 1:length(clasifNMF)){
    a[j]<-etiquetas_train[clasifNMF[j],3]
    b[j]<-etiquetas_train[clasifLNMF[j],3]
  }
  grafica4[1, contador]<-sum(a==as.numeric(as.factor(etiquetas_test[,3])))/length(a)
  grafica4[2, contador]<-sum(b==as.numeric(as.factor(etiquetas_test[,3])))/length(b)
  contador<-contador+1
}

plot(c(10,15,20,25,30,35,40,45),grafica3[1,], type="l", main="Test set", sub="Clasificación file", col=rgb(0.2,0.4,0.1,0.7) , lwd=3 , pch=17,xlab="Rango",ylab="Nivel de precisión",ylim=c(min(grafica3[1,],grafica3[2,]),max(grafica3[1,],grafica3[2,])))
lines(c(10,15,20,25,30,35,40,45),grafica3[2,],col=rgb(0.9,0.4,0.1,0.7) , lwd=3 , pch=19 , type="b")
legend("topleft", 
       legend = c("NMF", "Local NMF"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.9,0.4,0.1,0.7)), 
       pch = c(17,19,19), 
       bty = "n", 
       pt.cex = 1, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.001, 0.001))
table(a,etiquetas_test[,3])
table(b,etiquetas_test[,3])
