##You can run this file to get the FIG2 in Section 4.2.
##You need to load the data in the folder.
library(ggplot2)
#####Wheat Protein
load("~/Documents/research/Joint Model/weekly portfolio/final_draft3/realdata/Wheat/N10stdXY.RData")
Pro_10<-matrix(0,200,4)
Pro_10[,4]<-Outp[,2];Pro_10[,1]<-Outp[,5];Pro_10[,2]<-Outp[,7];Pro_10[,3]<-Outp[,6];
load("~/Documents/research/Joint Model/weekly portfolio/final_draft3/realdata/Wheat/N30stdXY.RData")
Pro_30<-matrix(0,200,4)
Pro_30[,4]<-Outp[,2];Pro_30[,1]<-Outp[,5];Pro_30[,2]<-Outp[,7];Pro_30[,3]<-Outp[,6];
PE=c(as.vector(Pro_10),as.vector(Pro_30))
Method=c(rep(c('Plug-In','IN','PACE','Proposed'),each=200),rep(c('Plug-In','IN','PACE','Proposed'),each=200))
N =c(rep('N=10',800), rep('N=30',800))
Data.df<-data.frame(Method,N,PE)
Data.df$Method <- factor(Data.df$Method , levels=c('Plug-In','IN','PACE','Proposed'))
ggplot(Data.df, aes(x=Method, y=PE, fill=N)) + ylab("Prediction error")+
  geom_boxplot() +
  scale_y_log10(breaks=c(0.1,0.2,0.4,0.6,0.8,1,10),limits=c(0.1,10))



