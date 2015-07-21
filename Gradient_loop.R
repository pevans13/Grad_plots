#script to loop through different models

#load up packages
library(ggplot2)
library(plyr)
library(lme4)
library(MuMIn)

Gradient<- read.csv("FinPlots7.csv")

pdf("Figures/Histogram.pdf")
names<-colnames(Gradient[,5:ncol(Gradient)])
for (i in 1:(ncol(Gradient)-4)){
  #this produces a histogram for each of your variables
  Hist<-ggplot(Gradient,aes(x=Gradient[[4+i]]))+geom_histogram()+xlab(names[i])+
    ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Hist)
}
dev.off()

for (i in 1:(ncol(Gradient)-4)){
  Shapiro_wilk<-do.call("rbind", with(Gradient6, tapply(Gradient[[4+2]], Plot,
  function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  Shapiro_wilk<-data.frame(Shapiro_wilk)
  if (length(Shapiro_wilk[Shapiro_wilk$p.value <= 0.05, ])>0){
    Modnull1<-glmer(Gradient[[4+2]]~ 1 +(1|Site),data=Gradient)
    Mod_cat<-glmer(Gradient[[4+2]] ~ Plot + (1| Site), data = Gradient)
    Mod_cont<-glmer(Gradient[[4+2]] ~ SBAPC + (1| Site), data = Gradient)
    Mod_cont_NL<-glmer(Gradient[[4+2]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gradient)
    Model_sel<-model.sel(Modnull1,Mod_cat,Mod_cont,extra = r.squaredGLMM)
    
  }
}
