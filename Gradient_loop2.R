#script to loop through different models

rm(list = ls())

#load up packages
library(ggplot2)
library(plyr)
library(lme4)
library(MuMIn)
library(reshape)

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

#load data
Gradient<- read.csv("FinPlots99.csv")
str(Gradient)
Names <- read.csv("F:/PhD/Chapter 1 Gradient Plots/R Scripts/Grad_plots/Names.csv", header=FALSE, dec=",")

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

# Create another data.frame just for the variables that need converting
GradT<-subset(Gradient, select=c(Site,Plot, SBAPC, Depth,Fagus.deadwood..m3., Lying.Deadwood.Total..m3.,NO3O,NH4O,NH4M, PMNO, PMNM,K,Ca,Mg,Na,Al,Mn,Fe,P,MF_Org,EC,Soil_Temp,MFSR,C_Storage..tonnes.hectares.))
GradSQ<-subset(Gradient, select=c(Site,Plot, SBAPC, Depth,Fagus.deadwood..m3., Lying.Deadwood.Total..m3.,NO3O,NH4O,NH4M, PMNO, PMNM,K,Ca,Mg,Na,Al,Mn,Fe,P,MF_Org,EC,Soil_Temp,MFSR,C_Storage..tonnes.hectares.))
names(GradT)
# Log
logcols<-c("Depth","NH4O","Al")
GradT$Mn<-log(GradT$Mn+1)
GradT$Fe<-log(GradT$Fe+1)
GradT$Al
GradT[logcols] <- log(GradT[logcols]+1)
# Square root
sqcols<-c("Fagus.deadwood..m3.","Lying.Deadwood.Total..m3.","NO3O","NH4O","NH4M","PMNO","PMNM","K","Ca","Na","MF_Org","Soil_Temp","MFSR")
GradT[sqcols] <- (GradT[sqcols]^0.5)
GradSQ<-subset(GradT,select=c(Site,SBAPC,Plot,Fagus.deadwood..m3.,Lying.Deadwood.Total..m3.,
                              NO3O,NH4O,NH4M,PMNO,PMNM,K,Ca,Na,MF_Org,Soil_Temp,MFSR)) # Save as a separate data.frame to use for models later
GradLog<-subset(GradT,select=c(Site,SBAPC,Plot,Depth,NH4O,Al,Mn,Fe)) # Save as a separate data.frame to use for models later

write.csv(GradLog,"Log variables for models.csv")
Lognames<-(Names[1:8,15])
for (i in 4:(ncol(GradLog))){
  Shapiro_wilk<-do.call("rbind", with(GradLog, tapply(GradLog[[i]], Plot,
                                                    function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Lognames[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:(ncol(GradLog))){
  b<-bartlett.test(resid(lm(GradLog[[i]]~Plot))~Plot,data=GradLog) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Lognames[i],".doc",sep="BT"))
}

for (i in 4:(ncol(GradLog)))
{
  Modnull1<-lmer(GradLog[[i]]~ 1 +(1|Site),data=GradLog)
  Mod_cont<-lmer(GradLog[[i]] ~ SBAPC + (1| Site), data = GradLog)
  Mod_cont_NL<-lmer(GradLog[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = GradLog)
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(Lognames[i],".csv",sep="Mod"),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=GradLog$Site)
  Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA)^2)
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=GradLog,aes(x=SBAPC*100,y=(GradLog[[i]]^2),colour=Site,group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred),groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(Lognames[i])+ggtitle(Lognames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}

# For capsule sites
GradCap<-subset(Gradient, select=c(Site,Plot, SBAPC, RCNit, RCAmm))
# Delete NA sites
GradCap<-completeFun(GradCap, "RCAmm")
sqcols<-c("RCAmm","RCNit")
GradCap[sqcols] <- sqrt(GradCap[sqcols])
str(GradCap)
# Inverse (1/x)
invcols<-c("Mg","EC","P","C_Storage..tonnes.hectares.")
GradT[invcols] <- 1/(GradT[invcols])

### Normality and hetero of transformed ariables, based on All_results.xls ##
names(GradT)
write.csv(GradT,"Transformed variables.csv")
Trannames<-(Names[1:24,8])
for (i in 4:(ncol(GradT))){
  Shapiro_wilk<-do.call("rbind", with(GradT, tapply(GradT[[i]], Plot,
                                                       function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Trannames[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:(ncol(GradT))){
  b<-bartlett.test(resid(lm(GradT[[i]]~Plot))~Plot,data=GradT) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Trannames[i],".doc",sep="BT"))
}



names(GradCap)
write.csv(GradCap,"Capsule variables.csv")
Capnames<-(Names[1:5,9])
for (i in 4:(ncol(GradCap))){
  Shapiro_wilk<-do.call("rbind", with(GradCap, tapply(GradCap[[i]], Plot,
                                                    function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Capnames[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:(ncol(GradCap))){
  b<-bartlett.test(resid(lm(GradCap[[i]]~Plot))~Plot,data=GradCap) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Capnames[i],".doc",sep="BT"))
}

pdf("Figures/Gradient_models_Square root variables.pdf")
## Use only gradSQ (variables that have been square-rooted)
write.csv(GradSQ,"Squared variables for models.csv")
SQnames<-(Names[1:16,11])
names(GradSQ)
for (i in 4:(ncol(GradSQ)))
{
  Modnull1<-lmer(GradSQ[[i]]~ 1 +(1|Site),data=GradSQ)
  Mod_cont<-lmer(GradSQ[[i]] ~ SBAPC + (1| Site), data = GradSQ)
  Mod_cont_NL<-lmer(GradSQ[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = GradSQ)
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(SQnames[i],".csv",sep="Mod"),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=GradSQ$Site)
  Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA)^2)
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=GradSQ,aes(x=SBAPC*100,y=(GradSQ[[i]]^2),colour=Site,group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred),groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(SQnames[i])+ggtitle(SQnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}
dev.off()


# Save as dataframe the important information that will always be used (i.e. site, sbapc, SBA, Plot and Soil_Type)
Grad<-(Gradient[c(1:5)])
# save a dataframe that includes the variables that will need transforming via qlogis
names(Gradient[c(102:128)])
Gradql<-Gradient[c(1:5,102:128)]
# convert them with qlogis
# CompleteFun before a variable signifies that NA existed in those variables and are therefore removed
# NA remove function

### LINE DOESN'T WORK - qlcols<-c("X100LOIO","X100LOIM","X100Clay","X100Sand","X100Silt","X100Bracken","X100Moss","X100BareGround","X100BG.M","X100Litter")
Gradql$X100LOIO<-qlogis(Gradql$X100LOIO)
Gradql$X100LOIM<-qlogis(Gradql$X100LOIM)
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
# Subset conditions which never have any values for the 5th stage (total collapse)
GradCond<-subset(Overall_Cond, select=c(Site,Plot, SBAPC, X100Crown_Loss,X100Leave_Loss,X100Crown_Condition,X100Discolouration,X100Condition_unhealthiness))
names(GradCond)
# Write the files to transfer the exact names to the Names csv
write.csv(Gradql,"Gradql.csv")
write.csv(GradCond,"GradCond.csv")
# Read in file with names to display on the y-axis

pdf("Figures/Gradient_models_Qlogis.pdf")
#loop to produce different plots for all the variables based on best model
QLnames<-(Names[1:23,5])
names(Gradql[1:23])
names(Gradql)
write.csv(Gradql,"Qlogis variables.csv")
for (i in 6:(ncol(Gradql))){
  Shapiro_wilk<-do.call("rbind", with(Gradql, tapply(Gradql[[i]], Plot,
                                                     function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(QLnames[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:(ncol(Gradql))){
  b<-bartlett.test(resid(lm(Gradql[[i]]~Plot))~Plot,data=Gradql) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(QLnames[i],".doc",sep="BT"))
}

for (i in 6:(ncol(Gradql)))
{
  Modnull1<-lmer(Gradql[[i]]~ 1 +(1|Site),data=Gradql)
  Mod_cont<-lmer(Gradql[[i]] ~ SBAPC + (1| Site), data = Gradql)
  Mod_cont_NL<-lmer(Gradql[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gradql)
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(QLnames[i],".csv",sep="Mod"),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradql$Site)
  Preds$PredR<-(predict(model.avg(Model_sel),newdata =Preds))
  Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=Gradql,aes(x=SBAPC*100,y=(plogis(Gradql[[i]]))*100,colour=Site,group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(plogis(Pred))*100,groups=NULL),colour="black",size=2) # Tried with plogis - line is always too high
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(QLnames[i])+ggtitle(QLnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}
dev.off()

# Condition of beech trees
pdf("Figures/Gradient_models_Qlogis_Cond.pdf")
Condnames<-(Names[1:8,6])
Condnames
head(GradCond)
for (i in 4:(ncol(GradCond)))
{
  Modnull1<-lmer(GradCond[[i]]~ 1 +(1|Site),data=GradCond)
  Mod_cont<-lmer(GradCond[[i]] ~ SBAPC + (1| Site), data = GradCond)
  Mod_cont_NL<-lmer(GradCond[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = GradCond)
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(Condnames[i],".csv",sep=""),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=GradCond$Site)
  Preds$PredR<-(predict(model.avg(Model_sel),newdata =Preds))
  Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=GradCond,aes(x=SBAPC*100,y=(plogis(GradCond[[i]]))*100,colour=Site,group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(plogis(Pred))*100,groups=NULL),colour="black",size=2) # Tried with plogis - line is always too high
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(Condnames[i])+ggtitle(Condnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}
dev.off()


pdf("Figures/Gradient_models_Poisson.pdf")
# set the names of the variables thast are going to be modelled with a poisson error distribution structure
names(Gradient[1:33])
Poivar <- Gradient[c(1:8,13:33)] # variables that need a poisson error distribution
Poivar <- subset(Poivar, select=-c(Fagusno,FagusJuv,FagusSL,QuerSL))
names(Poivar)
Pnames<-(Names[1:25,10])
write.csv(Poivar,"Poisson variables.csv")
for (i in 6:(ncol(Poivar))){
  Shapiro_wilk<-do.call("rbind", with(Poivar, tapply(Poivar[[i]], Plot,
                                                      function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Pnames[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:(ncol(Poivar))){
  b<-bartlett.test(resid(lm(Poivar[[i]]~Plot))~Plot,data=Poivar) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Pnames[i],".doc",sep="BT"))
}

for (i in 6:(ncol(Poivar))){
      Modnull1<-glmer(Poivar[[i]]~ 1 +(1|Site),data=Poivar,family="poisson")
      Mod_cont<-glmer(Poivar[[i]] ~ SBAPC + (1| Site), data = Poivar,family="poisson")
      Mod_cont_NL<-glmer(Poivar[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Poivar,family="poisson")
      Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
      print(Model_sel)
      write.csv(Model_sel,paste(Pnames[i],".csv",sep="Mod"),row.names=T)
      Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Poivar$Site)
      Preds$PredR<-exp(predict(model.avg(Model_sel),newdata =Preds))
      Preds$Pred<-exp(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
      theme_set(theme_bw(base_size=12))
      Plot1<-ggplot(data=Poivar,aes(x=SBAPC*100,y=Poivar[[i]],group=Site))+geom_point(size=2)
      Plot2<-Plot1+geom_line(data=Preds,aes(y=Pred,groups=NULL),colour="black",size=2)
      Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
      Plot4<-Plot3+xlab("Percentage loss of Basal area")+ylab(Pnames[i])+ggtitle(Pnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
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
Gusvar<-subset(Gusvar,select=-c(C_Storage..tonnes.hectares.,K,Ca,Mg,Na,Al,Mn,Fe,P,MF_Org,Fagus.deadwood..m3.,Lying.Deadwood.Total..m3.,
                                EC,RCAmm,RCNit,Soil_Temp,MFSR,Depth,LOIM,LOIO))
write.csv(Gusvar,"Gusvar.csv")
names(Gusvar)
Gusnames<-(Names[1:26,14])
for (i in 6:(ncol(Gusvar))){
  Modnull1<-lmer(Gusvar[[i]]~ 1 +(1|Site),data=Gusvar)
  Mod_cont<-lmer(Gusvar[[i]] ~ SBAPC + (1| Site), data = Gusvar)
  Mod_cont_NL<-lmer(Gusvar[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gusvar)
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(Gusnames[i],".csv",sep=""),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gusvar$Site)
  Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=Gusvar,aes(x=SBAPC*100,y=Gusvar[[i]],group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred),groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab(Gusnames[i])+ggtitle(Gusnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}
dev.off()

## The rest ##
Gradrest<-subset(Gradient, select=c(Site,Plot, SBAPC, Vsml,Med,Lrg,Huge,Fagusno,FagusJuv,QuerSL,
                                DBH,Height,SnagTot,SnagFagus,SnagIlex,SnagQuer))
Gradrest1<-subset(Gradient, select=c(Site,Plot, SBAPC, Vsml,Med,Lrg,Huge,Fagusno,DBH,Height,FagusSL))
Gradrest2<-completeFun(Gradrest1,"Height")

names(Gradrest)
write.csv(Gradrest,"The Rest!.csv")
write.csv(Gradrest2,"The Rest2.csv")
Restnames<-(Names[1:16,12]) # Names for the rest
for (i in 13:15){
  Shapiro_wilk<-do.call("rbind", with(Gradrest, tapply(Gradrest[[i]], Plot,
                                                     function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Restnames[i],".csv",sep="SW"),row.names=T)
}
for (i in 13:15){
  b<-bartlett.test(resid(lm(Gradrest[[i]]~Plot))~Plot,data=Gradrest) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Restnames[i],".doc",sep="BT"))
}

for (i in 4:16){
  Modnull1<-lmer(Gradrest[[i]]~ 1 +(1|Site),data=Gradrest)
  Mod_cont<-lmer(Gradrest[[i]] ~ SBAPC + (1| Site), data = Gradrest)
  Mod_cont_NL<-lmer(Gradrest[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gradrest)
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(Restnames[i],".csv",sep="Mod"),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradrest$Site)
  Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=Gradrest,aes(x=SBAPC*100,y=Gradrest[[i]],group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred),groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab(Restnames[i])+ggtitle(Restnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}

Restnames2<-(Names[1:10,13]) # Names for the rest
for (i in 4:(ncol(Gradrest2))){
  Shapiro_wilk<-do.call("rbind", with(Gradrest2, tapply(Gradrest2[[i]], Plot,
                                                       function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Restnames2[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:(ncol(Gradrest))){
  b<-bartlett.test(resid(lm(Gradrest2[[i]]~Plot))~Plot,data=Gradrest2) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Restnames2[i],".doc",sep="BT"))
}

for (i in 4:(ncol(Gradrest2))){
  Modnull1<-lmer(Gradrest2[[i]]~ 1 +(1|Site),data=Gradrest2)
  Mod_cont<-lmer(Gradrest2[[i]] ~ SBAPC + (1| Site), data = Gradrest2)
  Mod_cont_NL<-lmer(Gradrest2[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gradrest2)
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(Restnames2[i],".csv",sep="Mod"),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradrest2$Site)
  Preds$Pred<-exp(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=Gradrest2,aes(x=SBAPC*100,y=Gradrest2[[i]],group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=log(Pred),groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab(Restnames2[i])+ggtitle(Restnames2[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}


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
Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab("Soil respiration rate")+ggtitle("Soil respiration rate")+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
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

########### NO3M ##############
names(Gradient)
Modnull1<-lmer(NO3M~ 1 +(1|Site),data=Gradient)
Mod_cont<-lmer(NO3M ~ SBAPC +(1| Site), data = Gradient)
Mod_cont_NL<-lmer(NO3M~SBAPC +(1| Site)+I(SBAPC^2)+(1| Site), data = Gradient)
Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
print(Model_sel)
write.csv(Model_sel,"NO3M.csv",row.names=T)
Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradient$Site )
# Prodcue confidence intervals
Preds$NO3M <- predict(Mod_cont_NL,Preds,re.form=NA)
mm <- model.matrix(terms(Mod_cont_NL),Preds)
pvar1 <- diag(mm %*% tcrossprod(vcov(Mod_cont_NL),mm))
tvar1 <- pvar1+VarCorr(Mod_cont_NL)$Site[1]
cmult <- 2
Preds <- data.frame(
  Preds
  , plo = Preds$NO3M-cmult*sqrt(pvar1)
  , phi = Preds$NO3M+cmult*sqrt(pvar1)
  , tlo = Preds$NO3M-cmult*sqrt(tvar1)
  , thi = Preds$NO3M+cmult*sqrt(tvar1)
)

Preds$Pred<-exp(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
theme_set(theme_bw(base_size=12))
Plot1<-ggplot(data=Gradient,aes(x=SBAPC*100,y=NO3M,group=Site))+geom_point(size=2)
Plot2<-Plot1+geom_line(data=Preds,aes(y=log(Pred),groups=NULL),colour="black",size=2)+geom_ribbon(data=Preds,aes(ymax=phi,ymin=plo),alpha=0.01,colour=NA)
Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab("Nitrates in mineral layer")+ggtitle("Nitrates in mineral layer")+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
print(Plot4)

dev.off()

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
for (i in 6:(ncol(Poivar))){
theme_set(theme_bw(base_size=12))
m<-ggplot(subset(Gradient_summary, variable %in% c((Pnames[[i]])),aes(x=Plot,y=Mean))
m2<-m+geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.1,size=1.3)+geom_line(size=1)+geom_point(size=10,shape=20,col="black")+
  theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))+
  theme(axis.title.y = element_text(size = rel(2), angle = 90),
        axis.title.x = element_text(size = rel(2)))
m2
m3<-m+geom_pointrange(aes(ymin=Mean+SE,ymax=Mean-SE), width = 1)+geom_line(group=1)
m3
print(m2)
}
dev.off()

Fungi_s<-subset(Gradient_summary, variable==names[[i]])
Plot1<-ggplot(data=Fungi_s,aes(x=Plot,y=Mean),group=1)+geom_point()+geom_pointrange(aes(ymin=Mean+SE,ymax=Mean-SE), width = 1)+geom_line(group=1)


######## Poisson variables ######
names(Gradient[1:33])
Poivar <- Gradient[c(1:8,13:33)] # variables that need a poisson error distribution
names(Poivar)
Pnames<-(Names[1:33,2])

for (i in 6:(ncol(Poivar))){
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=Poivar,aes(x=Plot,y=Poivar[[i]])+geom_point(size=2))
 
}
dev.off()

############ other stats #############
### Explore gradient ##
names(Gradient)
Gradnames<-(Names[1:128,7])
names(Gradient[63:74])
for (i in 63:74){
  Shapiro_wilk<-do.call("rbind", with(Gradient, tapply(Gradient[[i]], Plot,
                                                       function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Gradnames[i],".csv",sep="SW"),row.names=T)
}

names(Gradient[30:45])
Shapiro_wilk<-do.call("rbind", with(Gradient, tapply(Lying.Deadwood.Total..m3., Plot,
                                                     function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
write.csv(Shapiro_wilk,"Biomass.csv",row.names=T)

for (i in 63:74){
  b<-bartlett.test(resid(lm(Gradient[[i]]~Plot))~Plot,data=Gradient) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Gradnames[i],".doc",sep="BT"))
}

Poivar <- Gradient[c(1:8,13:33)] # variables that need a poisson error distribution
names(Poivar)
for (i in 6:12){
  Shapiro_wilk<-do.call("rbind", with(Poivar, tapply(Poivar[[i]], Plot,
                                                    function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Pnames[i],".csv",sep="SW"),row.names=T)
}

for (i in 6:(ncol(Poivar))){
  b<-bartlett.test(resid(lm(Poivar[[i]]~Plot))~Plot,data=Poivar) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Pnames[i],".doc",sep="BT"))
}

# Trial to test loop
bartlett.test(resid(lm(Sward~Plot))~Plot,data=Poivar) 

for (i in 16:29){
  Shapiro_wilk<-do.call("rbind", with(Poivar, tapply(Poivar[[i]], Plot,
                                                     function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Pnames[i],".csv",sep="SW"),row.names=T)
}

# Trial a single to confirm loop worked
Shapiro_wilk<-do.call("rbind", with(Poivar, tapply(FagusSL, Plot,
                                                   function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
write.csv(Shapiro_wilk,paste("VVV.csv",sep="SW"),row.names=T)


Grad2<-Gradient[ , colSums(is.na(Gradient)) == 0] # Exclude columns (variables) that have NAs
write.csv(Grad2,"Grad2.csv")
names(Grad2)
namesSW<-Names[1:96,3]
for (i in 6:9){
  Shapiro_wilk<-do.call("rbind", with(Grad2, tapply(Grad2[[i]], Plot,
                                              function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(namesSW[i],".csv",sep="SW"),row.names=T)
}

for (i in 14){
  Shapiro_wilk<-do.call("rbind", with(Grad2, tapply(Grad2[[i]], Plot,
                                                    function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(namesSW[i],".csv",sep="SW"),row.names=T)
}

names(Grad2)
for (i in 17:27){
  Shapiro_wilk<-do.call("rbind", with(Grad2, tapply(Grad2[[i]], Plot,
                                                    function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(namesSW[i],".csv",sep="SW"),row.names=T)
}


### GradQL ##
QLnames
for (i in 6:(ncol(Gradql))){
  Shapiro_wilk<-do.call("rbind", with(Gradql, tapply(Gradql[[i]], Plot,
                                                    function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(QLnames[i],".csv",sep="SW"),row.names=T)
}

### GradCond ##
for (i in 4:(ncol(GradCond))){
  Shapiro_wilk<-do.call("rbind", with(GradCond, tapply(GradCond[[i]], Plot,
                                                     function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Condnames[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:(ncol(GradCond))){
  b<-bartlett.test(resid(lm(GradCond[[i]]~Plot))~Plot,data=GradCond) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Condnames[i],".doc",sep="BT"))
}

# Trial a single to confirm loop worked
Shapiro_wilk<-do.call("rbind", with(Gradient, tapply(log(Al), Plot,
                                                   function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
write.csv(Shapiro_wilk,paste("VVV.csv",sep="SW"),row.names=T)



for (i in 6:126){
ll<-summary(aov(Gradient[[i]]~Plot,data=Gradient))
capture.output(ll,file=paste(names[i],".doc",sep="AOV"))
}

for (i in 6:126){
pp<-TukeyHSD(aov(Gradient[[i]]~Plot,data=Gradient)) # Tukey post-hoc test 
capture.output(pp,file=paste(names[i],".doc",sep="Tuk"))
}