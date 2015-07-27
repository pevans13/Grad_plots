#script to loop through different models

rm(list = ls())

#load up packages
library(ggplot2)
library(plyr)
library(lme4)
library(MuMIn)
library(reshape)

se <- function(x) sqrt(var(x)/length(x))

#load data
Gradient<- read.csv("FinPlots7.csv")
str(Gradient)


#loop to produce different histograms for all the variables
pdf("Figures/Histogram.pdf")
names<-colnames(Gradient[,5:ncol(Gradient)])
for (i in 5:(ncol(Gradient))){
  #this produces a histogram for each of your variables
  Hist<-ggplot(Gradient,aes(x=Gradient[[4+i]]))+geom_histogram()+xlab(names[i])+
    ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Hist)
}
dev.off()

<<<<<<< HEAD
#loop to produce graphs and standard error bars
pdf("Figures/Standard error plot graphs.pdf")
names<-colnames(Gradient[,5:ncol(Gradient)])
for (i in 5:(ncol(Gradient))){
  #this produces a graph with standard error bars
  SE<-ggplot(Gradient,aes(x=Gradient[[4+i]]))+geom_line()+geom_errorbar()+xlab(names[i])+
    ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(SE)
}

for (i in 10:(ncol(Gradient))){
newSE <- summarySE(Gradient, measurevar=names[i], groupvars=c("Plot"))
}

g<-ggplot(Gradient, aes(x=Plot, y=names[i],group=1)) + 
  geom_errorbar(aes(ymin=names[i]-se, ymax=names[i]+se), width=0.1,size=1.3) +
  geom_line(size=1)+geom_point(size=10,shape=20,col="black")
g

  # Change the axis text
  g2<-g + theme(axis.text.x=element_text(angle=55, size=14, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=14, vjust=0.5))+
    labs(x="Stage of collapse", y=names[i])
  g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
  g4<-g3+theme(axis.title.y = element_text(size = rel(1), angle = 90,vjust=1.5),
               axis.title.x = element_text(size = rel(1)))
  g4
  # Change the aesthetics
  GF1<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
  GF1
}
dev.off()


=======
>>>>>>> 28b234acda44917f5a82a74bce14d808e9fcc5bc
pdf("Figures/Gradient_models.pdf")
#loop to go through different variables and produce models and figures for each of these
for (i in 5:(ncol(Gradient))){
    if(is.integer(Gradient[[i]])==T){
    Modnull1<-glmer(Gradient[[i]]~ 1 +(1|Site),data=Gradient,family="poisson")
    Mod_cont<-glmer(Gradient[[i]] ~ SBAPC + (1| Site), data = Gradient,family="poisson")
    Mod_cont_NL<-glmer(Gradient[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gradient,family="poisson")
    Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
    Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradient$Site)
    Preds$PredR<-exp(predict(model.avg(Model_sel),newdata =Preds))
    Preds$Pred<-exp(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
    theme_set(theme_bw(base_size=12))
    Plot1<-ggplot(data=Gradient,aes(x=SBAPC,y=Gradient[[i]],colour=Site,group=Site))+geom_point(size=4,shape=1)
    Plot2<-Plot1+geom_line(data=Preds,aes(y=PredR),size=0.5)+geom_line(data=Preds,aes(y=Pred,groups=NULL),colour="black",size=2)
    Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
    Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(names[i])+ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
    print(Plot4)
    }
 else{
   Modnull1<-glmer(Gradient[[i]]~ 1 +(1|Site),data=Gradient)
   Mod_cont<-glmer(Gradient[[i]] ~ SBAPC + (1| Site), data = Gradient)
   Mod_cont_NL<-glmer(Gradient[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gradient)
   Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
   Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradient$Site)
   Preds$PredR<-(predict(model.avg(Model_sel),newdata =Preds))
   Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
   theme_set(theme_bw(base_size=12))
   Plot1<-ggplot(data=Gradient,aes(x=SBAPC,y=Gradient[[i]],colour=Site,group=Site))+geom_point(size=4,shape=1)
   Plot2<-Plot1+geom_line(data=Preds,aes(y=PredR),size=0.5)+geom_line(data=Preds,aes(y=Pred,groups=NULL),colour="black",size=2)
   Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
   Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(names[i])+ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
   print(Plot4)
 }
}
dev.off()
<<<<<<< HEAD
=======

#loop to estimate standard errors of plots for each variable
for (i in 5:(ncol(Gradient))){
  print(i)
  summarySE(Gradient, measurevar=names[[i]], groupvars=c("Plot"))
}

head(Gradient)

Gradient_melt<-melt(Gradient,id.vars = c("Site","Plot"))
Gradient_melt$value<-ifelse(is.na(Gradient_melt$value),0,Gradient_melt$value)
Gradient_summary<-ddply(Gradient_melt,.(Plot,variable),summarise,Mean=mean(value),SE=se(value))



ggplot(Gradient_summary,aes(x=Plot,y=Mean,ymax=Mean+(2*SE),ymin=Mean-(2*SE)))+geom_pointrange()+geom_line(group=1)+facet_wrap(~variable,scales = "free_y")
ggsave("Figures/Gradient_summary.pdf",height=15,width=30,dpi=500,units="in")
head(Gradient_melt)
>>>>>>> 28b234acda44917f5a82a74bce14d808e9fcc5bc
