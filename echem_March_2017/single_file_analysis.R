library(foreign)
library(psych)
<<<<<<< HEAD

cleanFile<- function(filename){
=======
library(foreign)


cleanFile<- function(filename){
  steps <- matrix(data = NA, nrow =6, ncol=7)
  colnames(steps) <- c('name','AA','choline1','choline2','choline3','H2O2','slope')
>>>>>>> origin/master
  data1 <-read.csv(filename, col.names = c("marker","sent1","enz1","sent2","enz2"))
  data1 <- na.omit(data1)
  data1<-data.frame(data1)
  data1<-data1[-c(1:6),] #remove the first 6 rows, which are garbage
<<<<<<< HEAD
  data1[,2:5]<-data1[,2:5]*-1
  # data1[,2:5]<-data1[,2:5]*1000
  
  for(y in 2:5){
    H2O2in <- which(data1$marker=="5")
    H2O2done <- which(data1$marker=="6")
    H2O2Response <- data1[H2O2done,y]-data1[H2O2in,y]
    data1[y] <- data1[y]/H2O2Response}
  
  data1$electrode1 <- data1$enz1-data1$sent1
  data1$electrode2 <- data1$enz2-data1$sent2
  
  for(y in 2:7){
    name <- colnames(data1)[y]
    plot(data1[,y], type = "l", main = name)}
  data1 <- data.frame(data1)
  
  return(data1)
  }

makeSteps <- function(data1){
  steps <- matrix(data = NA, nrow =6, ncol=7)
  colnames(steps) <- c('name','AA','Ch1','Ch2','Ch3','H2O2','slope')
=======
  # data1[,2:5]<-data1[,2:5]*-1
  # data1[,2:5]<-data1[,2:5]*1000
>>>>>>> origin/master
  AAin <- which(data1$marker=="1")
  Ch1in <- which(data1$marker=="2")
  Ch2in <- which(data1$marker=="3")
  Ch3in <- which(data1$marker=="4")
  H2O2in <- which(data1$marker=="5")
  H2O2done <- which(data1$marker=="6")
<<<<<<< HEAD
=======
  
#   for(x in 2:5){
#     bsln <- median(data1[1:AAin,x])
#     data1[,x]<-data1[,x]-bsln
#   }
  
  data1$electrode1 <- data1$enz1-data1$sent1
  data1$electrode2 <- data1$enz2-data1$sent2
  
>>>>>>> origin/master
  z=1
  for(y in 2:7) {
    #creates matrix with steps for each channel
    name <- colnames(data1)[y]
    AA <- data1[Ch1in,y]-data1[AAin,y]
    Ch1 <- data1[Ch2in,y]-data1[Ch1in,y]
    Ch2 <- data1[Ch3in,y]-data1[Ch2in,y]
    Ch3 <- data1[H2O2in,y]-data1[Ch3in,y]
    H2O2 <- data1[H2O2done,y]-data1[H2O2in,y]
<<<<<<< HEAD
=======
    
    slopeMat1 <- matrix(data=NA,nrow = 4,ncol = 2)
    colnames(slopeMat1) <- c('uM','pA')
    slopeMat1[,1]<-c(0,20,40,60)
    slopeMat1[,2]<-c(data1[Ch1in,y],data1[Ch2in,y],data1[Ch3in,y],data1[H2O2in,y])
    slopeMat1 <- data.frame(slopeMat1)
    fit1 <- lm(pA ~ uM, data=slopeMat1)
    slope1 <- fit1$coefficients[[2]]      
    
    steps[z,]<-c(name,AA,Ch1,Ch2,Ch3,H2O2,slope1)
>>>>>>> origin/master
    
    current <-c(data1[Ch1in,y] , data1[Ch2in,y] , data1[Ch3in,y] , data1[H2O2in,y])
    choline <- c(0 , 20 , 40 , 60)
    fit <- lm(current~choline)
    plot(current~choline)
    mySlope <- fit$coefficients[[2]]
    
    steps[z,]<-c(name,AA,Ch1,Ch2,Ch3,H2O2, mySlope)
    
    z=z+1}
<<<<<<< HEAD
  
  steps <- data.frame(steps)
  steps$selectivity <- (125 * steps$slope)/steps$AA
  steps$AA <-  as.character(steps$AA)
  steps$AA <- as.numeric(steps$AA)
  steps$Ch1 <-  as.character(steps$Ch1)
  steps$Ch1 <- as.numeric(steps$Ch1)
  steps$Ch2 <-  as.character(steps$Ch2)
  steps$Ch2 <- as.numeric(steps$Ch2)
  steps$Ch3 <-  as.character(steps$Ch3)
  steps$Ch3 <- as.numeric(steps$Ch3)
  steps$H2O2 <-  as.character(steps$H2O2)
  steps$H2O2 <- as.numeric(steps$H2O2)
  steps$slope <- as.character(steps$slope)
  steps$slope <- as.numeric(steps$slope)
  steps$selectivity <- as.character(steps$selectivity)
  steps$selectivity <- as.numeric(steps$selectivity)
  
  
  return(steps)}

getSlope <- function(data1, electrodeColNum){
  AAin <- which(data1$marker=="1")
  Ch1in <- which(data1$marker=="2")
  Ch2in <- which(data1$marker=="3")
  Ch3in <- which(data1$marker=="4")
  H2O2in <- which(data1$marker=="5")
  H2O2done <- which(data1$marker=="6")
  
  current <-c(data1[Ch1in,electrodeColNum] , data1[Ch2in,electrodeColNum] , data1[Ch3in,electrodeColNum] , data1[H2O2in,electrodeColNum])
  choline <- c(0 , 20 , 40 , 60)
  fit <- lm(current~choline)
  plot(current~choline)
  mySlope <- fit$coefficients[[2]]
  return(mySlope) }



=======
  steps<-data.frame(steps)
  return(steps)
>>>>>>> origin/master


<<<<<<< HEAD
=======

mySlope = function(dataframe, column){
  data1<- dataframe
  slopeMat1 <- matrix(data=NA,nrow = 4,ncol = 2)
  colnames(slopeMat1) <- c('uM','pA')
  slopeMat1[,1]<-c(0,20,40,60)
  slopeMat1[,2]<-x#c(data1[Ch1in,column],data1[Ch2in,column],data1[Ch3in,column],data1[H2O2in,column])
  slopeMat1 <- data.frame(slopeMat1)
  fit1 <- lm(pA ~ uM, data=slopeMat1)
  slope1 <- fit1$coefficients[[2]]
  return(slope1)}

>>>>>>> origin/master

