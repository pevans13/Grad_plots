#script to loop through different models

rm(list = ls())

#load up packages
library(ggplot2)
library(plyr)
library(lme4)
library(MuMIn)
library(reshape)

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

#load data
Gradient<- read.csv("FinPlots99.csv")
str(Gradient)

# Delete the rows with no trees present
GradSt14 <- Gradient[-c(5, 10, 15,20,25,30,35,40,45,50,55,60), ] 
# Split the DBH sizes into equal quantile sizes
quantile(GradSt14$DBH,probs=seq(0,1,0.25))


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


pdf("Figures/Gradient_models.pdf")
#loop to go through different variables and produce models and figures for each of these
str(Gradient)
names(Gradient)

Gradient<- read.csv("FinPlots99.csv")
names(Gradient[c(1:7)])

# Save as dataframe the important information that will always be used (i.e. site, sbapc, SBA, Plot and Soil_Type)
Grad<-(Gradient[c(1:5)])
# save a dataframe that includes the variables that will need transforming via qlogis
names(Gradient[c(102:126)])
Gradql<-Gradient[c(1:5,102:126)]
# convert them with qlogis
# CompleteFun before a variable signifies that NA existed in those variables and are therefore removed
# NA remove function
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

Gradql$X100Clay<-qlogis(Gradql$X100Clay)
Gradql$X100Sand<-qlogis(Gradql$X100Sand)
Gradql$X100Silt<-qlogis(Gradql$X100Silt)
Gradql$X100Bracken<-qlogis(Gradql$X100Bracken)
Gradql$X100Moss<-qlogis(Gradql$X100Moss)
Gradql$X100BareGround<-qlogis(Gradql$X100BareGround)
Gradql$X100BG.M<-qlogis(Gradql$X100BG.M)
Gradql$X100Litter<-qlogis(Gradql$X100Litter)
Fagbro<-completeFun(Gradql,"X100FagusBrowse")
Fagbro$X100FagusBrowse<-qlogis(Fagbro$X100FagusBrowse) # dataframe saved as Fagbro
Ilbro<-completeFun(Gradql,"X100IlexBrowse")
Ilbro$X100IlexBrowse<-qlogis(Ilbro$X100IlexBrowse) # dataframe saved as Ilbro
Gradql$X100Trampling<-qlogis(Gradql$X100Trampling)
Gradql$X100Hollyshrub<-qlogis(Gradql$X100Hollyshrub)
Gradql$X100Rubusshrub<-qlogis(Gradql$X100Rubusshrub)
Gradql$X100IlexEat<-qlogis(Gradql$X100IlexEat)
Eatrub<-completeFun(Gradql,"X100RubEat")
Eatrub$X100RubEat<-qlogis(Eatrub$X100RubEat) # dataframe saved as Eatrub
Overall_Cond<-completeFun(Gradql,"X100Crown_Loss")
Overall_Cond$X100Crown_Loss<-qlogis(Overall_Cond$X100Crown_Loss)
# the dataframe saved here, Overall_Cond, can be used for most condition metrics as there is no beech presence in only stage 5
Overall_Cond$X100Leave_Loss<-qlogis(Overall_Cond$X100Leave_Loss)
Overall_Cond$X100Crown_Condition<-qlogis(Overall_Cond$X100Crown_Condition)
Overall_Cond$X100Discolouration<-qlogis(Overall_Cond$X100Discolouration)
Overall_Cond$X100Condition_unhealthiness<-qlogis(Overall_Cond$X100Condition_unhealthiness)
Understorey<-completeFun(Gradql,"X100Undercond")
Understorey$X100Undercond<-qlogis(Understorey$X100Undercond) # dataframe: Understorey
Gradql$X100Grass<-qlogis(Gradql$X100Grass)
Gradql$X100Canopy_openness<-qlogis(Gradql$X100Canopy_openness)
Gradql$X100Under_openness<-qlogis(Gradql$X100Under_openness)
Gradql$X100Canopy_Open_Total<-qlogis(Gradql$X100Canopy_Open_Total)
# delete the columns that are now in other dataframes
Gradql<-subset(Gradql, select=-c(X100FagusBrowse,X100IlexBrowse,X100RubEat,X100Crown_Loss,X100Leave_Loss,X100Crown_Condition,X100Discolouration,X100Condition_unhealthiness,X100Undercond))

# Read in file with names to display on the y-axis
Names <- read.csv("F:/PhD/Chapter 1 Gradient Plots/R Scripts/Grad_plots/Names.csv", header=FALSE, dec=",")
rownames(Names)

pdf("Figures/Gradient_models_Qlogis.pdf")
#loop to produce different plots for all the variables based on best model
names(Gradql)
## names<-colnames(Gradql[,1:21])
names<-(Names[97:126,2])
names
for (i in 6:(ncol(Gradql)))
{
  Modnull1<-lmer(Gradql[[i]]~ 1 +(1|Site),data=Gradql)
  Mod_cont<-lmer(Gradql[[i]] ~ SBAPC + (1| Site), data = Gradql)
  Mod_cont_NL<-lmer(Gradql[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gradql)
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(names[i],".csv",sep=""),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradql$Site)
  Preds$PredR<-exp(predict(model.avg(Model_sel),newdata =Preds))
  Preds$Pred<-exp(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=Gradql,aes(x=SBAPC*100,y=plogis(Gradql[[i]])*100,colour=Site,group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred)*100,groups=NULL),colour="black",size=2) # Tried with plogis - line is always too high
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(names[i])+ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}
dev.off()

pdf("Figures/Gradient_models_Poisson.pdf")
# set the names of the variables thast are going to be modelled with a poisson error distribution structure
names(Gradient[1:33])
Poivar <- Gradient[c(1:8,13:33)] # variables that need a poisson error distribution
names(Poivar)
names<-colnames(Poivar[,1:29])
names<-(Names[1:33,2])
for (i in 6:(ncol(Poivar))){
      Modnull1<-glmer(Poivar[[i]]~ 1 +(1|Site),data=Poivar,family="poisson")
      Mod_cont<-glmer(Poivar[[i]] ~ SBAPC + (1| Site), data = Poivar,family="poisson")
      Mod_cont_NL<-glmer(Poivar[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Poivar,family="poisson")
      Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
      print(Model_sel)
      write.csv(Model_sel,paste(names[i],".csv",sep=""),row.names=T)
      Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Poivar$Site)
      Preds$PredR<-exp(predict(model.avg(Model_sel),newdata =Preds))
      Preds$Pred<-exp(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
      theme_set(theme_bw(base_size=12))
      Plot1<-ggplot(data=Poivar,aes(x=SBAPC*100,y=Poivar[[i]],group=Site))+geom_point(size=2)
      Plot2<-Plot1+geom_line(data=Preds,aes(y=Pred,groups=NULL),colour="black",size=2)
      Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
      Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(names[i])+ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
      print(Plot4)
              }
dev.off()

# Run models where herbivory may have an effect - Sward height and seedlings
pdf("Figures/Gradient_models_Poisson_Dung.pdf")
names(Gradient[1:33])
Pdung <- Gradient[c(1:6,24,30:33)] # variables that need a poisson error distribution
names(Pdung)
names<-colnames(Pdung[,1:11])
for (i in 8:(ncol(Pdung))){
  Modnull1<-glmer(Pdung[[i]]~ 1 +(1|Site),data=Pdung,family="poisson")
  Mod_cont<-glmer(Pdung[[i]] ~ SBAPC:DungTot + (1| Site), data = Pdung,family="poisson")
  Mod_cont_NL<-glmer(Pdung[[i]]~SBAPC:DungTot+I(SBAPC^2)+(1| Site), data = Pdung,family="poisson")
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(names[i],".csv",sep=""),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Pdung$Site,DungTot=mean(Pdung$DungTot))
  Preds$PredR<-exp(predict(model.avg(Model_sel),newdata =Preds))
  Preds$Pred<-exp(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=Pdung,aes(x=SBAPC*100,y=Pdung[[i]],group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=Pred,groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(names[i])+ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}
dev.off()

# Now the remaining Gaussian plots
pdf("Figures/Gradient_models_Gaussian.pdf")
# set the names of the variables 
names(Gradient)
Gusvar <- Gradient[c(1:5,34:74)]
names(Gusvar)
head(Gusvar)
names<-colnames(Gusvar[1:46])
for (i in 6:(ncol(Gusvar))){
  Modnull1<-lmer(Gusvar[[i]]~ 1 +(1|Site),data=Gusvar)
  Mod_cont<-lmer(Gusvar[[i]] ~ SBAPC + (1| Site), data = Gusvar)
  Mod_cont_NL<-lmer(Gusvar[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gusvar)
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(names[i],".csv",sep=""),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gusvar$Site)
  Preds$Pred<-exp(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=Gusvar,aes(x=SBAPC*100,y=Gusvar[[i]],group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=log(Pred),groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(names[i])+ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}
dev.off()

# Soil Respiration separately
Modnull1<-lmer(Soil_Respiration~ 1 +(1|Site),data=Gradient)
Mod_cont<-lmer(Soil_Respiration ~ SBAPC +Soil_Temp+Amb_Temp+MFSR +(1| Site), data = Gradient)
Mod_cont_NL<-lmer(Soil_Respiration~SBAPC +Soil_Temp+Amb_Temp+MFSR +(1| Site)+I(SBAPC^2)+(1| Site), data = Gradient)
Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
print(Model_sel)
write.csv(Model_sel,"SR.csv",row.names=T)
Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradient$Site,Soil_Temp=mean(Gradient$Soil_Temp),Amb_Temp=mean(Gradient$Amb_Temp),MFSR=mean(Gradient$MFSR) )
# Prodcue confidence intervals
Preds$Soil_Respiration <- predict(Mod_cont_NL,Preds,re.form=NA)
mm <- model.matrix(terms(Mod_cont_NL),Preds)
pvar1 <- diag(mm %*% tcrossprod(vcov(Mod_cont_NL),mm))
tvar1 <- pvar1+VarCorr(Mod_cont_NL)$Site[1]
cmult <- 2
Preds <- data.frame(
  Preds
  , plo = Preds$Soil_Respiration-cmult*sqrt(pvar1)
  , phi = Preds$Soil_Respiration+cmult*sqrt(pvar1)
  , tlo = Preds$Soil_Respiration-cmult*sqrt(tvar1)
  , thi = Preds$Soil_Respiration+cmult*sqrt(tvar1)
)

Preds$Pred<-exp(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
theme_set(theme_bw(base_size=12))
Plot1<-ggplot(data=Gradient,aes(x=SBAPC*100,y=Soil_Respiration,group=Site))+geom_point(size=2)
Plot2<-Plot1+geom_line(data=Preds,aes(y=log(Pred),groups=NULL),colour="black",size=2)+geom_ribbon(data=Preds,aes(ymax=phi,ymin=plo),alpha=0.01,colour=NA)
Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab("Soil respiration rate")+ggtitle("Soil respiration rate")+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
print(Plot4)

########### Biomass ##############
names(Gradient)
Modnull1<-lmer(Biomass..m3.ha.~ 1 +(1|Site),data=Gradient)
Mod_cont<-lmer(Biomass..m3.ha. ~ SBAPC +(1| Site), data = Gradient)
Mod_cont_NL<-lmer(Biomass..m3.ha.~SBAPC +(1| Site)+I(SBAPC^2)+(1| Site), data = Gradient)
Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
print(Model_sel)
write.csv(Model_sel,"SR.csv",row.names=T)
Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradient$Site )
# Prodcue confidence intervals
Preds$Biomass..m3.ha. <- predict(Mod_cont_NL,Preds,re.form=NA)
mm <- model.matrix(terms(Mod_cont_NL),Preds)
pvar1 <- diag(mm %*% tcrossprod(vcov(Mod_cont_NL),mm))
tvar1 <- pvar1+VarCorr(Mod_cont_NL)$Site[1]
cmult <- 2
Preds <- data.frame(
  Preds
  , plo = Preds$Biomass..m3.ha.-cmult*sqrt(pvar1)
  , phi = Preds$Biomass..m3.ha.+cmult*sqrt(pvar1)
  , tlo = Preds$Biomass..m3.ha.-cmult*sqrt(tvar1)
  , thi = Preds$Biomass..m3.ha.+cmult*sqrt(tvar1)
)

Preds$Pred<-exp(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
theme_set(theme_bw(base_size=12))
Plot1<-ggplot(data=Gradient,aes(x=SBAPC*100,y=Biomass..m3.ha.,group=Site))+geom_point(size=2)
Plot2<-Plot1+geom_line(data=Preds,aes(y=log(Pred),groups=NULL),colour="black",size=2)+geom_ribbon(data=Preds,aes(ymax=phi,ymin=plo),alpha=0.01,colour=NA)
Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab("Above ground biomass (m3 ha-1)")+ggtitle("Above ground biomass (m3 ha-1)")+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
print(Plot4)

############### Standard error ###############
#loop to estimate standard errors of plots for each variable
for (i in 5:(ncol(Gradient))){
  print(i)
  summarySE(Gradient, measurevar=names[[i]], groupvars=c("Plot"))
}

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
se <- function(x) sqrt(var(x)/length(x))
head(Gradient)

Gradient_melt<-melt(Gradient,id.vars = c("Site","Plot"))

Gradient_summary<-ddply(Gradient_melt,.(Plot,variable),summarise,Mean=mean(value,na.rm = T),SE=se(value))
head(Gradient_summary)

ggplot(Gradient_summary,aes(x=Plot,y=Mean,ymax=Mean+(2*SE),ymin=Mean-(2*SE)))+geom_pointrange()+geom_line(group=1)+facet_wrap(~variable,scales = "free_y")
ggsave("Figures/Gradient_summary.pdf",height=15,width=30,dpi=500,units="in")
head(Gradient_melt)

theme_set(theme_bw(base_size=12))
Fungi_s<-subset(Gradient_summary, variable=='Fungi')
Plot1<-ggplot(data=Fungi_s,aes(x=Plot,y=Mean),group=1)+geom_point()+geom_pointrange(aes(ymin=Mean+SE,ymax=Mean-SE), width = 1)+geom_line(group=1)


pdf("Figures/Anova_Poisson.pdf")
# set the names of the variables thast are going to be modelled with a standard errors
theme_set(theme_bw(base_size=12))
m<-ggplot(subset(Gradient_summary, variable %in% c("Sward")),aes(x=Plot,y=mean))
m<-m+geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.1,size=1.3)+geom_line(size=1)+geom_point(size=10,shape=20,col="black")
m
m<-m+geom_pointrange(aes(ymin=Mean+SE,ymax=Mean-SE), width = 1)+geom_line(group=1)
m

Fungi_s<-subset(Gradient_summary, variable==names[[i]])
Plot1<-ggplot(data=Fungi_s,aes(x=Plot,y=Mean),group=1)+geom_point()+geom_pointrange(aes(ymin=Mean+SE,ymax=Mean-SE), width = 1)+geom_line(group=1)

names(Gradient[1:33])
Poivar <- Gradient[c(1:8,13:33)] # variables that need a poisson error distribution
names(Poivar)
names<-(Names[1:33,2])

for (i in 6:(ncol(Poivar))){
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=Poivar,aes(x=Plot,y=Poivar[[i]])+geom_point(size=2))
 
}
dev.off()

############ other stats #############
names<-colnames(Gradient[,1:29])
for (i in 6:10){
  Shapiro_wilk<-do.call("rbind", with(Gradient, tapply(names[i], Plot,
                                                       function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(names[i],".csv",sep=""),row.names=T)
}

for (i in 6:10){
ll<-summary(aov(names[[i]]~Plot,data=Gradient)) # run if two tests are above 0.05. See one.way below, if not
}

by(Gradient,Gradient$Plot,function(x){
  anova(lm(names[[i]]~Plot,data=x))
})