#script to loop through different models

#load up packages
library(ggplot2)
library(plyr)
library(lme4)
library(MuMIn)

#load data
Gradient<- read.csv("FinPlots7.csv")
str(Gradient)


#loop to produce different histograms for all the variables
pdf("Figures/Histogram.pdf")
names<-colnames(Gradient[,5:ncol(Gradient)])
for (i in 1:(ncol(Gradient)-4)){
  #this produces a histogram for each of your variables
  Hist<-ggplot(Gradient,aes(x=Gradient[[4+i]]))+geom_histogram()+xlab(names[i])+
    ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Hist)
}
dev.off()


pdf("Figures/Gradient_models.pdf")
#loop to go through different variables and produce models and figures for each of these
for (i in 1:(ncol(Gradient)-4)){
    if(is.integer(Gradient[[i]])==T){
    Modnull1<-glmer(Gradient[[i]]~ 1 +(1|Site),data=Gradient,family="poisson")
    Mod_cont<-glmer(Gradient[[i]] ~ SBAPC + (1| Site), data = Gradient,family="poisson")
    Mod_cont_NL<-glmer(Gradient[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gradient,family="poisson")
    Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
    Preds<-data.frame(SBAPC=Gradient$SBAPC,
                      Site=Gradient$Site,
            PredR=exp(predict(model.avg(Model_sel))),
            Pred=exp(predict(model.avg(Model_sel),re.form=NA)))
    
    theme_set(theme_bw(base_size=12))
    Plot1<-ggplot(data=Gradient,aes(x=SBAPC,y=Gradient[[i]],colour=Site,group=Site))+geom_point(size=2)
    Plot2<-Plot1+geom_line(data=Preds,aes(y=PredR),size=0.5)+geom_line(data=Preds,aes(y=Pred,groups=NULL),colour="black",size=2)
    Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
    Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(names[i])+ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
    }
 else{
   Modnull1<-glmer(Gradient[[i]]~ 1 +(1|Site),data=Gradient)
   Mod_cont<-glmer(Gradient[[i]] ~ SBAPC + (1| Site), data = Gradient)
   Mod_cont_NL<-glmer(Gradient[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gradient)
   Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
   Preds<-data.frame(SBAPC=Gradient$SBAPC,
                     Site=Gradient$Site,
                     PredR=exp(predict(model.avg(Model_sel))),
                     Pred=exp(predict(model.avg(Model_sel),re.form=NA)))
   theme_set(theme_bw(base_size=12))
   Plot1<-ggplot(data=Gradient,aes(x=SBAPC,y=Gradient[[i]],colour=Site,group=Site))+geom_point(size=2)
   Plot2<-Plot1+geom_line(data=Preds,aes(y=PredR),size=0.5)+geom_line(data=Preds,aes(y=Pred,groups=NULL),colour="black",size=2)
   Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
   Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(names[i])+ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
 }
}
dev.off()