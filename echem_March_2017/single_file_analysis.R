library(foreign)
setwd("E:/echem_March_2017/CSVs")
directory <-list.files("E:/echem_March_2017/CSVs")

print(directory)

steps <- matrix(data = NA, nrow =6, ncol=6)
colnames(steps) <- c('name','AA','choline1','choline2','choline3','H2O2')


cleanFile<- function(i){
  filename=directory[i]
  data1 <-read.csv(filename, col.names = c("marker","sent1","enz1","sent2","enz2"))
  data1 <- na.omit(data1)
  data1<-data.frame(data1)
  data1<-data1[-c(1:6),] #remove the first 6 rows, which are garbage
  data1[,2:5]<-data1[,2:5]*-1
  data1[,2:5]<-data1[,2:5]*1000
  AAin <- which(data1$marker=="1")
  Ch1in <- which(data1$marker=="2")
  Ch2in <- which(data1$marker=="3")
  Ch3in <- which(data1$marker=="4")
  H2O2in <- which(data1$marker=="5")
  H2O2done <- which(data1$marker=="6")
  
  for(x in 2:5){
    bsln <- median(data1[1:AAin,x])
    data1[,x]<-data1[,x]-bsln
  }
  
  data1$electrode1 <- data1$enz1-data1$sent1
  data1$electrode2 <- data1$enz2-data1$sent2
  
  z=1
  for(y in 2:7) {
    #creates matrix with steps for each channel
    name <- colnames(data1)[y]
    AA <- data1[Ch1in,y]-data1[AAin,y]
    Ch1 <- data1[Ch2in,y]-data1[Ch1in,y]
    Ch2 <- data1[Ch3in,y]-data1[Ch2in,y]
    Ch3 <- data1[H2O2in,y]-data1[Ch3in,y]
    H2O2 <- data1[H2O2done,y]-data1[H2O2in,y]
    steps[z,]<-c(name,AA,Ch1,Ch2,Ch3,H2O2)
    
    #plots each channel
    plot(data1[,y], type = "l", main = name)
    z=z+1}
  
  return(data1)

}

#analysis of step values
library(psych)
steps<-data.frame(steps)
describe(steps)
