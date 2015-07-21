#script to loop through different models

#load up packages
library(ggplot2)

Gradient<- read.csv("FinPlots7.csv")

names<-colnames(Gradient[,5:ncol(Gradient)])
for (1:(ncol(Gradient)-4)){
  Hist<-ggplot(Gradient,aes(x=Gradient[[4+1]]))+geom_histogram()+xlab(names[1])+
    ggtitle(names[1])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
}



