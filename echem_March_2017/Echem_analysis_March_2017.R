library(foreign)
setwd("E:/echem_March_2017/CSVs")
directory <-list.files("E:/echem_March_2017/CSVs")
L <- length(directory)

steps <- matrix(data = NA, nrow =L*2, ncol=7)
colnames(steps) <- c('name','AA','choline1','choline2','choline3','H2O2','slope_pA_uM')

i=1
i2=1
for (item in directory){
  i3=i2+1
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
    bsln <- mean(data1[,x])
    data1[,x]<-data1[,x]-bsln
  }
  
  data1$electrode1 <- data1$enz1-data1$sent1
  data1$electrode2 <- data1$enz2-data1$sent2
  
  cholineAdditions <- data1[1:H2O2in,]
  title=paste(filename, " / ",i2)
  plot(cholineAdditions$electrode1,col=4, type="l", main = title)
  lines(cholineAdditions$enz1, col=3)
  lines(cholineAdditions$sent1, col=2)
  
  title=paste(filename, " / ",i3)
  plot(cholineAdditions$electrode2,col=4, type="l", main = title)
  lines(cholineAdditions$enz2, col=3)
  lines(cholineAdditions$sent2, col=2)
  
  electrode1AA <- data1[Ch1in,6]-data1[AAin,6]
  electrode1Ch1 <- data1[Ch2in,6]-data1[Ch1in,6]
  electrode1Ch2 <- data1[Ch3in,6]-data1[Ch2in,6]
  electrode1Ch3 <- data1[H2O2in,6]-data1[Ch3in,6]
  electrode1H2O2 <- data1[H2O2done,6]-data1[H2O2in,6]
  
  slopeMat1 <- matrix(data=NA,nrow = 4,ncol = 2)
  colnames(slopeMat1) <- c('uM','pA')
  slopeMat1[,1]<-c(0,20,40,60)
  slopeMat1[,2]<-c(data1[Ch1in,6],data1[Ch2in,6],data1[Ch3in,6],data1[H2O2in,6])
  slopeMat1 <- data.frame(slopeMat1)
  fit1 <- lm(pA ~ uM, data=slopeMat1)
  slope1 <- fit1$coefficients[[2]]
  
  steps[i2,]<-c(filename,electrode1AA,electrode1Ch1,electrode1Ch2,electrode1Ch3,electrode1H2O2,slope1)
  
  electrode2AA <- data1[Ch1in,7]-data1[AAin,7]
  electrode2Ch1 <- data1[Ch2in,7]-data1[Ch1in,7]
  electrode2Ch2 <- data1[Ch3in,7]-data1[Ch2in,7]
  electrode2Ch3 <- data1[H2O2in,7]-data1[Ch3in,7]
  electrode2H2O2 <- data1[H2O2done,7]-data1[H2O2in,7]
  
  slopeMat2 <- matrix(data=NA,nrow = 4,ncol = 2)
  colnames(slopeMat2) <- c('uM','pA')
  slopeMat2[,1]<-c(0,20,40,60)
  slopeMat2[,2]<-c(data1[Ch1in,6],data1[Ch2in,6],data1[Ch3in,6],data1[H2O2in,6])
  slopeMat2 <- data.frame(slopeMat2)
  fit2 <- lm(pA ~ uM, data=slopeMat2)
  slope2 <- fit2$coefficients[[2]]
  
  steps[i3,] <- c(filename, electrode2AA, electrode2Ch1, electrode2Ch2,electrode2Ch3,electrode2H2O2,slope2)
  
  i<-i+1
  i2<-i2+2
}

#analysis of step values
library(psych)
steps<-data.frame(steps)
steps$name <- as.character(steps$name)
steps$AA <-as.character(steps$AA)
steps$choline1 <-as.character(steps$choline1)
steps$choline2 <-as.character(steps$choline2)
steps$choline3 <-as.character(steps$choline3)
steps$H2O2 <-as.character(steps$H2O2)
steps$AA <-as.numeric(steps$AA)
steps$choline1 <-as.numeric(steps$choline1)
steps$choline2 <-as.numeric(steps$choline2)
steps$choline3 <-as.numeric(steps$choline3)
steps$H2O2 <-as.numeric(steps$H2O2)
steps$slope_pA_uM <- as.character(steps$slope_pA_uM)
steps$slope_pA_uM <- as.numeric(steps$slope_pA_uM)
describe(steps)

#Here I browsed through the plots of each electrode's activity and identified those that
#   either (1) showed only noise, (2) had no evidence of choline steps
cleanedSteps <- steps#enter the row numbers of electrodes you want to remove
cleanedSteps$cholineAvrg <- rowMeans(cleanedSteps[,(3:5)])
cleanedSteps$selectivity <- ((cleanedSteps$cholineAvrg) / (cleanedSteps$AA))*(127.5/20)
cleanedSteps$selectivity2 <- ((cleanedSteps$cholineAvrg) / (abs(cleanedSteps$AA)))*(127.5/20)
#I'm measuring selectivity by looking at the ratio of average choline response to AA response, 
#   multiplied by the ration of AA concentration to choline concentration per addition
describe(cleanedSteps)


bestElectrodes <- subset(cleanedSteps, cleanedSteps$selectivity > 50)
describe(bestElectrodes)
PercentMakeCut <- 100*(nrow(cleanedSteps)/nrow(steps))
PercentMakeCut #this is the percentage of my electrodes that are suitable for use at the end


