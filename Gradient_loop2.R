#script to loop through different models

rm(list = ls())

#load up packages
library(ggplot2)
library(plyr)
library(lme4)
library(MuMIn)
library(reshape)
library(segmented)
library(SiZer)
library(ez)


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

#load data
Gradient<- read.csv("FinPlots99.csv")

## Create SBA mean graph
newSE <- summarySE(Gradient, measurevar="SBA", groupvars=c("Plot"))
g<-ggplot(newSE, aes(x=Plot, y=SBA,group=1)) + 
  geom_errorbar(aes(ymin=SBA-se, ymax=SBA+se), width=0.1,size=1.3) +
  geom_line(size=1)+geom_point(size=10,shape=20,col="black")
# Change the axis text
g2<-g + theme(axis.text.x=element_text(angle=55, size=20, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=20, vjust=0.5))+
  labs(x="Stage of collapse", y="Stand basal area (SBA) (m2 ha-1)")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = 30, angle = 90,vjust=1.5),
             axis.title.x = element_text(size = 30))
g4
# Change the aesthetics
GF1<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
GF1

str(Gradient)
Names <- read.csv("F:/PhD/Chapter 1 Gradient Plots/R Scripts/Grad_plots/Names.csv", header=FALSE, dec=",")

# Split the DBH sizes into equal quantile sizes
Size_quantiles<-quantile(Gradient$DBH,probs=seq(0,1,0.25),na.rm=T); Size_quantiles

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
invcols<-c("P","EC","C_Storage..tonnes.hectares.","Mg")
GradT[invcols] <- 1/(GradT[invcols])
GradSQ<-subset(GradT,select=c(Site,SBAPC,Plot,Fagus.deadwood..m3.,Lying.Deadwood.Total..m3.,
                              NO3O,NH4M,PMNO,PMNM,K,Ca,Na,MF_Org,Soil_Temp,MFSR)) # Save as a separate data.frame to use for models later
GradLog<-subset(GradT,select=c(Site,SBAPC,Plot,Depth,NH4O,Al,Mn,Fe)) # Save as a separate data.frame to use for models later
GradInv<-subset(GradT, select=c(Site,Plot, SBAPC, P,EC,C_Storage..tonnes.hectares.,Mg))
write.csv(GradLog,"Log variables for models.csv")
write.csv(GradInv,"Inverse variables for models.csv")
Lognames<-(Names[1:8,15])
for (i in 4:(ncol(GradLog))){
  Shapiro_wilk<-do.call("rbind", with(GradLog, tapply(GradLog[[i]], Plot,
                                                    function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Lognames[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:(ncol(GradLog))){
  b<-bartlett.test(resid(lm(GradLog[[i]]~Plot))~Plot,data=GradLog) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Lognames[i],".doc",sep="BT"))
  kw<-kruskal.test(GradLog[[i]]~Plot,data=GradLog)
  capture.output(kw,file=paste(Lognames[i],".doc",sep="KW"))
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
  Preds$Pred<-exp(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=GradLog,aes(x=SBAPC*100,y=exp(GradLog[[i]]) ,group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred),groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab(Lognames[i])+ggtitle(Lognames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}

Invnames<-(Names[1:7,16])
for (i in 4:(ncol(GradInv))){
  Shapiro_wilk<-do.call("rbind", with(GradInv, tapply(GradInv[[i]], Plot,
                                                      function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Invnames[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:(ncol(GradInv))){
  b<-bartlett.test(resid(lm(GradInv[[i]]~Plot))~Plot,data=GradInv) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Invnames[i],".doc",sep="BT"))
  kw<-kruskal.test(GradInv[[i]]~Plot,data=GradInv)
  capture.output(kw,file=paste(Invnames[i],".doc",sep="KW"))
}

# pdf("Figures/Gradient_models_inv.pdf")
for (i in 4:(ncol(GradInv)))
{
  Modnull1<-lmer(GradInv[[i]]~ 1 +(1|Site),data=GradInv)
  Mod_cont<-lmer(GradInv[[i]] ~ SBAPC + (1| Site), data = GradInv)
  Mod_cont_NL<-lmer(GradInv[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = GradInv)
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(Invnames[i],".csv",sep="Mod"),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=GradInv$Site)
  Preds$Pred<-1/(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=GradInv,aes(x=SBAPC*100,y=1/(GradInv[[i]]),group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred),groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab(Invnames[i])+ggtitle(Invnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}
# dev.off()

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

### Normality and hetero of transformed variables, based on All_results.xls ##
# pdf("Figures/Gradient_models_Min.pdf")
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
  kw<-kruskal.test(GradCap[[i]]~Plot,data=GradCap)
  capture.output(kw,file=paste(Capnames[i],".doc",sep="KW"))
}

for (i in 4:(ncol(GradCap)))
{
  Modnull1<-lmer(GradCap[[i]]~ 1 +(1|Site),data=GradCap)
  Mod_cont<-lmer(GradCap[[i]] ~ SBAPC + (1| Site), data = GradCap)
  Mod_cont_NL<-lmer(GradCap[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = GradCap)
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(Capnames[i],".csv",sep="Mod"),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=GradCap$Site)
  Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA)^2)
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=GradCap,aes(x=SBAPC*100,y=(GradCap[[i]]^2),group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred),groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab(Capnames[i])+ggtitle(Capnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}

# dev.off()
## Threshold analysis of capsule nutrients
GradCap<-subset(Gradient, select=c(Site,Plot, SBAPC, RCNit, RCAmm))
# Delete NA sites
GradCap<-completeFun(GradCap, "RCAmm");Capnames<-(Names[1:5,9])
for (i in 4:5){
  Thres1 <- piecewise.linear(GradCap$SBAPC, GradCap[[i]], middle=1, CI=TRUE,bootstrap.samples = 1000, sig.level = 0.05)
  capture.output(Thres1,file=paste(Capnames[i],".doc",sep="Thres"))
}

# pdf("Figures/Gradient_models_Square root variables.pdf")
## Use only gradSQ (variables that have been square-rooted)
write.csv(GradSQ,"Squared variables for models.csv")
SQnames<-(Names[1:15,11])
names(GradSQ)
for (i in 4:(ncol(GradSQ))){
  Shapiro_wilk<-do.call("rbind", with(GradSQ, tapply(GradSQ[[i]], Plot,
                                                      function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(SQnames[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:(ncol(GradSQ))){
  b<-bartlett.test(resid(lm(GradSQ[[i]]~Plot))~Plot,data=GradSQ) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(SQnames[i],".doc",sep="BT"))
  kw<-kruskal.test(GradSQ[[i]]~Plot,data=GradSQ)
  capture.output(kw,file=paste(SQnames[i],".doc",sep="KW"))
}

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
  Plot1<-ggplot(data=GradSQ,aes(x=SBAPC*100,y=(GradSQ[[i]]^2) ,group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred),groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab(SQnames[i])+ggtitle(SQnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}
# dev.off()

# save a dataframe that includes the variables that will need transforming via qlogis
names(Gradient[c(102:128)])
Gradql<-Gradient[c(1:5,102:128)]
# convert them with qlogis

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

Gradql$X100Trampling<-qlogis(Gradql$X100Trampling)
Gradql$X100Hollyshrub<-qlogis(Gradql$X100Hollyshrub)
Gradql$X100Rubusshrub<-qlogis(Gradql$X100Rubusshrub)
Gradql$X100IlexEat<-qlogis(Gradql$X100IlexEat)

Overall_Cond<-completeFun(Gradql,"X100Crown_Loss")
Overall_Cond$X100Crown_Loss<-qlogis(Overall_Cond$X100Crown_Loss)
# the dataframe saved here, Overall_Cond, can be used for most condition metrics as there is no beech presence in only stage 5
Overall_Cond$X100Leave_Loss<-qlogis(Overall_Cond$X100Leave_Loss)
Overall_Cond$X100Crown_Condition<-qlogis(Overall_Cond$X100Crown_Condition)
Overall_Cond$X100Discolouration<-qlogis(Overall_Cond$X100Discolouration)
Overall_Cond$X100Condition_unhealthiness<-qlogis(Overall_Cond$X100Condition_unhealthiness)

Gradql$X100Grass<-qlogis(Gradql$X100Grass)
Gradql$X100Canopy_openness<-qlogis(Gradql$X100Canopy_openness)
Gradql$X100Under_openness<-qlogis(Gradql$X100Under_openness)
Gradql$X100Canopy_Open_Total<-qlogis(Gradql$X100Canopy_Open_Total)

# Separate the variables with NAs in different locations
Fagbro<-completeFun(Gradql,"X100FagusBrowse")
Fagbro$X100FagusBrowse<-qlogis(Fagbro$X100FagusBrowse) # dataframe saved as Fagbro
Ilbro<-completeFun(Gradql,"X100IlexBrowse")
Ilbro$X100IlexBrowse<-qlogis(Ilbro$X100IlexBrowse) # dataframe saved as Ilbro
Eatrub<-completeFun(Gradql,"X100RubEat")
Eatrub$X100RubEat<-qlogis(Eatrub$X100RubEat) # dataframe saved as Eatrub
Understorey<-completeFun(Gradql,"X100Undercond")
Understorey$X100Undercond<-qlogis(Understorey$X100Undercond) # dataframe: Understorey
Shapiro_wilk<-do.call("rbind", with(Understorey, tapply(X100Undercond, Plot,function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
write.csv(Shapiro_wilk,paste("UnderstoreySW.csv",sep="SW"),row.names=T)
bartlett.test(resid(lm(X100Undercond~Plot))~Plot,data=Understorey) # Homogeneity of Variance of residuals

kw<-kruskal.test(X100Undercond~Plot,data=Understorey)
capture.output(kw,file="X100UndercondKW.doc")


# delete the columns that are now in other dataframes
Gradql<-subset(Gradql, select=-c(X100FagusBrowse,X100IlexBrowse,X100RubEat,X100Crown_Loss,X100Leave_Loss,X100Crown_Condition,X100Discolouration,X100Condition_unhealthiness,X100Undercond))

# Subset conditions which never have any values for the 5th stage (total collapse)
GradCond<-subset(Overall_Cond, select=c(Site,Plot, SBAPC, X100Crown_Loss,X100Leave_Loss,X100Crown_Condition,X100Discolouration,X100Condition_unhealthiness))
names(GradCond)

# Write the files to transfer the exact names to the Names csv
write.csv(Gradql,"Gradql.csv")
write.csv(GradCond,"GradCond.csv")
# Read in file with names to display on the y-axis

# pdf("Figures/Gradient_models_Qlogis.pdf")
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
  kw<-kruskal.test(Gradql[[i]]~Plot,data=Gradql)
  capture.output(kw,file=paste(QLnames[i],".doc",sep="KW"))
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
  Plot1<-ggplot(data=Gradql,aes(x=SBAPC*100,y=(plogis(Gradql[[i]]))*100 ,group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(plogis(Pred))*100,groups=NULL),colour="black",size=2) # Tried with plogis - line is always too high
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab(QLnames[i])+ggtitle(QLnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}
# dev.off()

# Condition of beech trees
# pdf("Figures/Gradient_models_Qlogis_Cond.pdf")
Condnames<-(Names[1:8,6])
Condnames

### GradCond ##
for (i in 4:(ncol(GradCond))){
  Shapiro_wilk<-do.call("rbind", with(GradCond, tapply(GradCond[[i]], Plot,
                                                       function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Condnames[i],".csv",sep="SW"),row.names=T)
  b<-bartlett.test(resid(lm(GradCond[[i]]~Plot))~Plot,data=GradCond) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Condnames[i],".doc",sep="BT"))
}

head(GradCond)
for (i in 4:(ncol(GradCond))){
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
  Plot1<-ggplot(data=GradCond,aes(x=SBAPC*100,y=(plogis(GradCond[[i]]))*100 ,group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(plogis(Pred))*100,groups=NULL),colour="black",size=2) # Tried with plogis - line is always too high
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab(Condnames[i])+ggtitle(Condnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}

for (i in 4:(ncol(GradCond))){
  b<-bartlett.test(resid(lm(GradCond[[i]]~Plot))~Plot,data=GradCond) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Condnames[i],".doc",sep="BT"))
  kw<-kruskal.test(GradCond[[i]]~Plot,data=GradCond)
  capture.output(kw,file=paste(Condnames[i],".doc",sep="KW"))
}
# dev.off()


# pdf("Figures/Gradient_models_Poisson.pdf")
# set the names of the variables thast are going to be modelled with a poisson error distribution structure
names(Gradient[1:33])
Poivar <- Gradient[c(1:8,13:33)] # variables that need a poisson error distribution
Poivar <- subset(Poivar, select=-c(Fagusno,FagusJuv,FagusSL,QuerSL,IlexSL,TotSL))
names(Poivar)
Pnames<-(Names[1:23,10])
write.csv(Poivar,"Poisson variables.csv")
for (i in 6:(ncol(Poivar))){
  Shapiro_wilk<-do.call("rbind", with(Poivar, tapply(Poivar[[i]], Plot,
                                                      function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Pnames[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:(ncol(Poivar))){
  b<-bartlett.test(resid(lm(Poivar[[i]]~Plot))~Plot,data=Poivar) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Pnames[i],".doc",sep="BT"))
  kw<-kruskal.test(Poivar[[i]]~Plot,data=Poivar)
  capture.output(kw,file=paste(Pnames[i],".doc",sep="KW"))
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
      Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab(Pnames[i])+ggtitle(Pnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
      print(Plot4)
              }



# dev.off()

# Run models where herbivory may have an effect - Sward height and seedlings
# pdf("Figures/Gradient_models_Poisson_Dung.pdf")
names(Gradient[1:33])
Pdung <- Gradient[c(1:6,24,30:33)] # variables that need a poisson error distribution
names(Pdung)
write.csv(Pdung,"Pdung.csv")
SLnames<-(Names[1:11,18])
for (i in 8:(ncol(Pdung))){
  Modnull1<-glmer(Pdung[[i]]~ 1 +(1|Site),data=Pdung,family="poisson")
  Mod_cont<-glmer(Pdung[[i]] ~ SBAPC:DungTot + (1| Site), data = Pdung,family="poisson")
  Mod_cont_NL<-glmer(Pdung[[i]]~SBAPC:DungTot+I(SBAPC^2)+(1| Site), data = Pdung,family="poisson")
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(SLnames[i],".csv",sep="Mod"),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Pdung$Site,DungTot=mean(Pdung$DungTot))
  Preds$PredR<-exp(predict(model.avg(Model_sel),newdata =Preds))
  Preds$Pred<-exp(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=Pdung,aes(x=SBAPC*100,y=Pdung[[i]],group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=Pred,groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab(SLnames[i])+ggtitle(SLnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}

for (i in 4:(ncol(Pdung))){
  b<-bartlett.test(resid(lm(Pdung[[i]]~Plot))~Plot,data=Pdung) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(SLnames[i],".doc",sep="BT"))
  kw<-kruskal.test(Pdung[[i]]~Plot,data=Pdung)
  capture.output(kw,file=paste(SLnames[i],".doc",sep="KW"))
}
# dev.off()

# Now the remaining Gaussian plots
# pdf("Figures/Gradient_models_Gaussian.pdf")
# set the names of the variables 
names(Gradient)
Gusvar <- Gradient[c(1:5,34:74)]
Gusvar<-subset(Gusvar,select=-c(C_Storage..tonnes.hectares.,K,Ca,Mg,Na,Al,Mn,Fe,P,MF_Org,Fagus.deadwood..m3.,Lying.Deadwood.Total..m3.,
                                EC,RCAmm,RCNit,Soil_Temp,MFSR,Depth,LOIM,LOIO,PMNM,PMNO,NO3O,SnagIlex,SnagTot,SnagQuer,SnagFagus,
                                DBH,Height,NH4O,NH4M,Soil_Respiration,Biomass..m3.ha.))
write.csv(Gusvar,"Gusvar.csv")
names(Gusvar)
Gusnames<-(Names[1:13,14])

for (i in 6:(ncol(Gusvar))){
  Shapiro_wilk<-do.call("rbind", with(Gusvar, tapply(Gusvar[[i]], Plot,
                                                     function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Gusnames[i],".csv",sep="SW"),row.names=T)
}
for (i in 6:(ncol(Gusvar))){
  b<-bartlett.test(resid(lm(Gusvar[[i]]~Plot))~Plot,data=Gusvar) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Gusnames[i],".doc",sep="BT"))
  kw<-kruskal.test(Gusvar[[i]]~Plot,data=Gusvar)
  capture.output(kw,file=paste(Gusnames[i],".doc",sep="KW"))
}

for (i in 6:(ncol(Gusvar))){
  Modnull1<-lmer(Gusvar[[i]]~ 1 +(1|Site),data=Gusvar)
  Mod_cont<-lmer(Gusvar[[i]] ~ SBAPC + (1| Site), data = Gusvar)
  Mod_cont_NL<-lmer(Gusvar[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gusvar)
  Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
  print(Model_sel)
  write.csv(Model_sel,paste(Gusnames[i],".csv",sep="Mod"),row.names=T)
  Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gusvar$Site)
  Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
  theme_set(theme_bw(base_size=12))
  Plot1<-ggplot(data=Gusvar,aes(x=SBAPC*100,y=Gusvar[[i]],group=Site))+geom_point(size=2)
  Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred),groups=NULL),colour="black",size=2)
  Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
  Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab(Gusnames[i])+ggtitle(Gusnames[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
  print(Plot4)
}
# dev.off()

## The rest ##
Gradrest<-subset(Gradient, select=c(Site,Plot, SBAPC, Vsml,Med,Lrg,Huge,Fagusno,FagusJuv,QuerSL,
                                DBH,Height,SnagTot,SnagFagus,SnagIlex,SnagQuer))
Gradrest1<-subset(Gradient, select=c(Site,Plot, SBAPC, Vsml,Med,Lrg,Huge,Fagusno,DBH,Height))
Gradrest2<-completeFun(Gradrest1,"Height")
Gradrest<-subset(Gradrest, select=c(Site,Plot, SBAPC,FagusJuv,QuerSL,
                                    SnagTot,SnagFagus,SnagIlex,SnagQuer)) # Once the NA have been removed can save again as gradrest
names(Gradrest)
names(Gradrest2)
write.csv(Gradrest,"The Rest!.csv")
write.csv(Gradrest2,"The Rest2.csv")
Restnames<-(Names[1:9,12]) # Names for the rest
for (i in 6:8){
  Shapiro_wilk<-do.call("rbind", with(Gradrest, tapply(Gradrest[[i]], Plot,
                                                     function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Restnames[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:9){
  b<-bartlett.test(resid(lm(Gradrest[[i]]~Plot))~Plot,data=Gradrest) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Restnames[i],".doc",sep="BT"))
  kw<-kruskal.test(Gradrest[[i]]~Plot,data=Gradrest)
  capture.output(kw,file=paste(Restnames[i],".doc",sep="KW"))
}

for (i in 4:9){
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

Restnames2<-(Names[1:11,13]) # Names for the rest
for (i in 4:(ncol(Gradrest2))){
  Shapiro_wilk<-do.call("rbind", with(Gradrest2, tapply(Gradrest2[[i]], Plot,
                                                       function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(Restnames2[i],".csv",sep="SW"),row.names=T)
}
for (i in 4:(ncol(Gradrest2))){
  b<-bartlett.test(resid(lm(Gradrest2[[i]]~Plot))~Plot,data=Gradrest2) # Homogeneity of Variance of residuals
  capture.output(b,file=paste(Restnames2[i],".doc",sep="BT"))
  kw<-kruskal.test(Gradrest2[[i]]~Plot,data=Gradrest2)
  capture.output(kw,file=paste(Restnames2[i],".doc",sep="KW"))
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

pdf("Figures/Gradient_models_ci.pdf")
# Soil Respiration separately
Shapiro_wilk<-do.call("rbind", with(Gradient, tapply(Soil_Respiration, Plot,
                                                      function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
write.csv(Shapiro_wilk,paste("SRSW.csv"),row.names=T)
b<-bartlett.test(resid(lm(Soil_Respiration~Plot))~Plot,data=Gradient) # Homogeneity of Variance of residuals
capture.output(b,file=paste("srBT.doc"))
kw<-TukeyHSD(aov(Soil_Respiration~Plot,data=Gradient))
capture.output(kw,file="Soil_RespirationPH.doc")

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

names(Gradient)
Modnull1<-lmer(CNRatioSoil~ 1 +(1|Site),data=Gradient)
Mod_cont<-lmer(CNRatioSoil ~ SBAPC +(1| Site), data = Gradient)
Mod_cont_NL<-lmer(CNRatioSoil~SBAPC  +(1| Site)+I(SBAPC^2)+(1| Site), data = Gradient)
Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
print(Model_sel)
write.csv(Model_sel,"SR.csv",row.names=T)
Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradient$Site )
# Prodcue confidence intervals
Preds$CNRatioSoil <- predict(Mod_cont_NL,Preds,re.form=NA)
mm <- model.matrix(terms(Mod_cont_NL),Preds)
pvar1 <- diag(mm %*% tcrossprod(vcov(Mod_cont_NL),mm))
tvar1 <- pvar1+VarCorr(Mod_cont_NL)$Site[1]
cmult <- 2
Preds <- data.frame(
  Preds
  , plo = Preds$CNRatioSoil-cmult*sqrt(pvar1)
  , phi = Preds$CNRatioSoil+cmult*sqrt(pvar1)
  , tlo = Preds$CNRatioSoil-cmult*sqrt(tvar1)
  , thi = Preds$CNRatioSoil+cmult*sqrt(tvar1)
)

Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
theme_set(theme_bw(base_size=12))
Plot1<-ggplot(data=Gradient,aes(x=SBAPC*100,y=CNRatioSoil,group=Site))+geom_point(size=2)
Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred),groups=NULL),colour="black",size=2)+geom_ribbon(data=Preds,aes(ymax=phi,ymin=plo),alpha=0.01,colour=NA)
Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab("Carbon nitrogen ratio")+ggtitle("Carbon nitrogen ratio")+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
print(Plot4)

names(Gradient)
Modnull1<-lmer(LichenRichBeech~ 1 +(1|Site),data=Gradient)
Mod_cont<-lmer(LichenRichBeech ~ SBAPC+ (1| Site), data = Gradient)
Mod_cont_NL<-lmer(LichenRichBeech~SBAPC +(1| Site)+I(SBAPC^2)+(1| Site), data = Gradient)
Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
print(Model_sel)
write.csv(Model_sel,"SR.csv",row.names=T)
Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradient$Site )
# Prodcue confidence intervals
Preds$LichenRichBeech <- predict(Mod_cont_NL,Preds,re.form=NA)
mm <- model.matrix(terms(Mod_cont_NL),Preds)
pvar1 <- diag(mm %*% tcrossprod(vcov(Mod_cont_NL),mm))
tvar1 <- pvar1+VarCorr(Mod_cont_NL)$Site[1]
cmult <- 2
Preds <- data.frame(
  Preds
  , plo = Preds$LichenRichBeech-cmult*sqrt(pvar1)
  , phi = Preds$LichenRichBeech+cmult*sqrt(pvar1)
  , tlo = Preds$LichenRichBeech-cmult*sqrt(tvar1)
  , thi = Preds$LichenRichBeech+cmult*sqrt(tvar1)
)

Preds$Pred<-predict(model.avg(Model_sel),newdata =Preds,re.form=NA)
theme_set(theme_bw(base_size=12))
Plot1<-ggplot(data=Gradient,aes(x=SBAPC*100,y=LichenRichBeech,group=Site))+geom_point(size=2)
Plot2<-Plot1+geom_line(data=Preds,aes(y=Pred,groups=NULL),colour="black",size=2)+geom_ribbon(data=Preds,aes(ymax=phi,ymin=plo),alpha=0.01,colour=NA)
Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab("Lichen richness on beech")+ggtitle("Lichen richness on beech")+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
print(Plot4)

names(Gradient)
Modnull1<-lmer(Fungi~ 1 +(1|Site),data=Gradient)
Mod_cont<-lmer(Fungi ~ SBAPC+ (1| Site), data = Gradient)
Mod_cont_NL<-lmer(Fungi~SBAPC +(1| Site)+I(SBAPC^2)+(1| Site), data = Gradient)
Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
print(Model_sel)
write.csv(Model_sel,"SR.csv",row.names=T)
Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradient$Site )
# Prodcue confidence intervals
Preds$Fungi <- predict(Mod_cont_NL,Preds,re.form=NA)
mm <- model.matrix(terms(Mod_cont_NL),Preds)
pvar1 <- diag(mm %*% tcrossprod(vcov(Mod_cont_NL),mm))
tvar1 <- pvar1+VarCorr(Mod_cont_NL)$Site[1]
cmult <- 2
Preds <- data.frame(
  Preds
  , plo = Preds$Fungi-cmult*sqrt(pvar1)
  , phi = Preds$Fungi+cmult*sqrt(pvar1)
  , tlo = Preds$Fungi-cmult*sqrt(tvar1)
  , thi = Preds$Fungi+cmult*sqrt(tvar1)
)

Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
theme_set(theme_bw(base_size=12))
Plot1<-ggplot(data=Gradient,aes(x=SBAPC*100,y=Fungi,group=Site))+geom_point(size=2)
Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred),groups=NULL),colour="black",size=2)+geom_ribbon(data=Preds,aes(ymax=phi,ymin=plo),alpha=0.01,colour=NA)
Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab("Mycorrhizal fungal richness")+ggtitle("Mycorrhizal fungal richness")+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
print(Plot4)

names(Gradient)
Modnull1<-lmer(C_Storage..tonnes.hectares.~ 1 +(1|Site),data=Gradient)
Mod_cont<-lmer(C_Storage..tonnes.hectares. ~ SBAPC+ (1| Site), data = Gradient)
Mod_cont_NL<-lmer(C_Storage..tonnes.hectares.~SBAPC +(1| Site)+I(SBAPC^2)+(1| Site), data = Gradient)
Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
print(Model_sel)
write.csv(Model_sel,"SR.csv",row.names=T)
Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradient$Site )
# Prodcue confidence intervals
Preds$C_Storage..tonnes.hectares. <- predict(Mod_cont_NL,Preds,re.form=NA)
mm <- model.matrix(terms(Mod_cont_NL),Preds)
pvar1 <- diag(mm %*% tcrossprod(vcov(Mod_cont_NL),mm))
tvar1 <- pvar1+VarCorr(Mod_cont_NL)$Site[1]
cmult <- 2
Preds <- data.frame(
  Preds
  , plo = Preds$C_Storage..tonnes.hectares.-cmult*sqrt(pvar1)
  , phi = Preds$C_Storage..tonnes.hectares.+cmult*sqrt(pvar1)
  , tlo = Preds$C_Storage..tonnes.hectares.-cmult*sqrt(tvar1)
  , thi = Preds$C_Storage..tonnes.hectares.+cmult*sqrt(tvar1)
)

Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
theme_set(theme_bw(base_size=12))
Plot1<-ggplot(data=Gradient,aes(x=SBAPC*100,y=C_Storage..tonnes.hectares.,group=Site))+geom_point(size=2)
Plot2<-Plot1+geom_line(data=Preds,aes(y=(Pred),groups=NULL),colour="black",size=2)+geom_ribbon(data=Preds,aes(ymax=phi,ymin=plo),alpha=0.01,colour=NA)
Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab("Carbon storage (tonnes ha-1)")+ggtitle("Carbon storage (tonnes ha-1)")+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
print(Plot4)


########### Biomass ##############
kw<-kruskal.test(Biomass..m3.ha.~Plot,data=Gradient)
capture.output(kw,file="Biomass..m3.ha.KW.doc")

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

########### Organic.C..tonnes.hectare. ##############
Shapiro_wilk<-do.call("rbind", with(Gradient, tapply(Organic.C..tonnes.hectare., Plot,
                                                     function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
write.csv(Shapiro_wilk,paste("SoilcarbonSW.csv"),row.names=T)
b<-bartlett.test(resid(lm(Organic.C..tonnes.hectare.~Plot))~Plot,data=Gradient) # Homogeneity of Variance of residuals
capture.output(b,file=paste("SoilcarbonBT.doc"))
kw<-kruskal.test(Organic.C..tonnes.hectare.~Plot,data=Gradient)
capture.output(kw,file="Organic.C..tonnes.hectare.KW.doc")

jj<-oneway.test(Organic.C..tonnes.hectare.~Plot,data=Gradient)
capture.output(jj,file="SoilcarbonWelanova.doc")

names(Gradient)
Modnull1<-lmer(Organic.C..tonnes.hectare.~ 1 +(1|Site),data=Gradient)
Mod_cont<-lmer(Organic.C..tonnes.hectare. ~ SBAPC +(1| Site), data = Gradient)
Mod_cont_NL<-lmer(Organic.C..tonnes.hectare.~SBAPC +(1| Site)+I(SBAPC^2)+(1| Site), data = Gradient)
Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
print(Model_sel)
write.csv(Model_sel,"Soilcarbon.csv",row.names=T)
Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradient$Site )
# Prodcue confidence intervals
Preds$Organic.C..tonnes.hectare. <- predict(Mod_cont_NL,Preds,re.form=NA)
mm <- model.matrix(terms(Mod_cont_NL),Preds)
pvar1 <- diag(mm %*% tcrossprod(vcov(Mod_cont_NL),mm))
tvar1 <- pvar1+VarCorr(Mod_cont_NL)$Site[1]
cmult <- 2
Preds <- data.frame(
  Preds
  , plo = Preds$Organic.C..tonnes.hectare.-cmult*sqrt(pvar1)
  , phi = Preds$Organic.C..tonnes.hectare.+cmult*sqrt(pvar1)
  , tlo = Preds$Organic.C..tonnes.hectare.-cmult*sqrt(tvar1)
  , thi = Preds$Organic.C..tonnes.hectare.+cmult*sqrt(tvar1)
)

Preds$Pred<-predict(model.avg(Model_sel),newdata =Preds,re.form=NA)
theme_set(theme_bw(base_size=12))
Plot1<-ggplot(data=Gradient,aes(x=SBAPC*100,y=Organic.C..tonnes.hectare.,group=Site))+geom_point(size=2)
Plot2<-Plot1+geom_line(data=Preds,aes(y=Pred,groups=NULL),colour="black",size=2)+geom_ribbon(data=Preds,aes(ymax=phi,ymin=plo),alpha=0.01,colour=NA)
Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Plot4<-Plot3+xlab("Percentage loss of basal area")+ylab("Soil organic carbon (tonnes ha-1)")+ggtitle("Soil organic carbon (tonnes ha-1)")+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
print(Plot4)
dev.off()

Thres1 <- piecewise.linear(Gradient$SBAPC,Gradient$PMNM, middle=1, CI=TRUE,bootstrap.samples = 1000, sig.level = 0.05)
capture.output(Thres1,file="PNMNM.doc")

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


# pdf("Figures/Anova_Poisson.pdf")
# set the names of the variables thast are going to be modelled with a standard errors
for (i in 6:(ncol(Poivar))){
theme_set(theme_bw(base_size=12))
<<<<<<< HEAD
m<-ggplot(subset(Gradient_summary, variable %in% c("Sward")),aes(x=Plot,y=Mean))+geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.1,size=1.3)+geom_line(size=1)+geom_point(size=10,shape=20,col="black")
m
m<-m+geom_pointrange(aes(ymin=Mean+SE,ymax=Mean-SE), width = 1)+geom_line(group=1)
m
=======
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
>>>>>>> 2790830b094413d0acd08c099bfcce3a897526c7

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
<<<<<<< HEAD
names2<-colnames(Gradient[,1:29])
for (i in 6:10){
  Shapiro_wilk<-do.call("rbind", with(Gradient, tapply(names[i], Plot,
=======
### Explore gradient ##
names(Gradient)
Gradnames<-(Names[1:128,7])
names(Gradient[63:74])
for (i in 63:74){
  Shapiro_wilk<-do.call("rbind", with(Gradient, tapply(Gradient[[i]], Plot,
>>>>>>> 2790830b094413d0acd08c099bfcce3a897526c7
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

<<<<<<< HEAD
i<-6


for (i in 6:10){
ll<-summary(aov(names2[[6]]~Plot,data=Gradient)) # run if two tests are above 0.05. See one.way below, if not
=======
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

### GradQL ##
QLnames
for (i in 6:(ncol(Gradql))){
  Shapiro_wilk<-do.call("rbind", with(Gradql, tapply(Gradql[[i]], Plot,
                                                    function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
  write.csv(Shapiro_wilk,paste(QLnames[i],".csv",sep="SW"),row.names=T)
>>>>>>> 2790830b094413d0acd08c099bfcce3a897526c7
}

# Trial a single to confirm loop worked
Shapiro_wilk<-do.call("rbind", with(Gradient, tapply(sqrt(PMNO), Plot,
                                                   function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
write.csv(Shapiro_wilk,paste("VVV.csv",sep="SW"),row.names=T)

for (i in 6:126){
pp<-TukeyHSD(aov(Gradient[[i]]~Plot,data=Gradient)) # Tukey post-hoc test 
capture.output(pp,file=paste(names[i],".doc",sep="Tuk"))
}

names(Gradient)
Threshold<-subset(Gradient,select=-c(Vsml,Med,Lrg,Huge,RCAmm,RCNit,DBH,Height,Squirrel,Debarking))
write.csv(Threshold,"TP.csv")
namesThres<-Names[1:91,17]
# for (i in 6:91){
# for (i in 25:91){
for (i in 88:91){
Thres1 <- piecewise.linear(Threshold$SBAPC, Threshold[[i]], middle=1, CI=TRUE,bootstrap.samples = 1000, sig.level = 0.05)
capture.output(Thres1,file=paste(namesThres[i],".doc",sep="Thres"))
}
Gradql<-Gradient[c(1:5,102:129)]
Overall_Cond<-completeFun(Gradql,"X100Crown_Loss"); Condnames<-(Names[1:8,6])
for (i in 4:81){
  Thres1 <- piecewise.linear(Overall_Cond$SBAPC, Overall_Cond[[i]], middle=1, CI=TRUE,bootstrap.samples = 1000, sig.level = 0.05)
  capture.output(Thres1,file=paste(Condnames[i],".doc",sep="Thres"))
}

## Anova tests based on previous S-W and bartlett test results ##
GradAnova<-subset(Gradient,select=c(Site,Plot,Fungi,LichenRicHol,Depth,NO3O,NH4O,PMNM,CNRatioSoil, BD..g.cm3.,K,Ca,Mg,Na,P,MF_Org,pH,EC,RCAmm,Soil_Respiration,
                                    X100Clay,X100Sand,X100Canopy_Open_Total,PMNO))
logcols<-c("Depth","NH4O")
invcol<-c("Mg","P","EC")
SQcols<-c("NO3O","PMNM","K","Ca","Na","MF_Org","RCAmm","PMNO")
Alogit<-c("X100Clay","X100Sand","X100Canopy_Open_Total")
GradAnova[logcols] <- log(GradAnova[logcols]+1)
GradAnova[SQcols] <- sqrt(GradAnova[SQcols])
GradAnova[invcol] <- 1/(GradAnova[invcol])
GradAnova[Alogit] <- (qlogis(GradAnova[Alogit]))
write.csv(GradAnova,"Anovavar.csv")
ANnames<-(Names[1:24,19])

for (i in 3:24){
  jj<-oneway.test(GradAnova[[i]]~Plot,data=GradAnova)
  capture.output(jj,file=paste(ANnames[i],".doc",sep="ANOVA"))
}

# Welch's correction anova
WelAnova<-subset(Gradient,select=c(Site,Plot,CompTot,CompGF,LichenRichBeech,Fagus.deadwood..m3.,Lying.Deadwood.Total..m3.,Al,
                                   Amb_Temp, C_Storage..tonnes.hectares.,X100Leave_Loss,X100Discolouration,X100Canopy_openness))
SQcols<-c("Fagus.deadwood..m3.","Lying.Deadwood.Total..m3.")
logcols<-c("Al")
invcol<-c("C_Storage..tonnes.hectares.")
WelAnova[logcols] <- log(WelAnova[logcols]+1)
WelAnova[SQcols] <- sqrt(WelAnova[SQcols])
WelAnova[invcol] <- 1/(WelAnova[invcol])
write.csv(WelAnova,"WelAnova.csv")
Welnames<-(Names[1:13,20])
for (i in 3:13){
  jj<-oneway.test(WelAnova[[i]]~Plot,data=WelAnova)
  capture.output(jj,file=paste(Welnames[i],".doc",sep="ANOVA"))
}
##### Games-Howell ####
# Carries out post-hoc tests for Welch's Anova
tukey <- function(  data,      	
                    group,					
                    method=c("Tukey", "Games-Howell"))	
{
  OK <- complete.cases(data, group)			
  data <- data[OK]
  group <- factor(group[OK])
  n <- tapply(data, group, length)		
  a <- length(n)					
  phi.e <- sum(n)-a					
  Mean <- tapply(data, group, mean)			
  Variance <- tapply(data, group, var)		
  result1 <- cbind(n, Mean, Variance)		
  rownames(result1) <- paste("Group", 1:a, sep="")
  method <- match.arg(method)
  if (method == "Tukey") {				# Tukey 
    v.e <- sum((n-1)*Variance)/phi.e		
    t <- combn(a, 2, function(ij)		
      abs(diff(Mean[ij]))/sqrt(v.e*sum(1/n[ij])) )
    p <- ptukey(t*sqrt(2), a, phi.e, lower.tail=FALSE)	
    Tukey <- cbind(t, p)					
    rownames(Tukey) <- combn(a, 2, paste, collapse=":")
    return(list(result1=result1, Tukey=Tukey, phi=phi.e, v=v.e))
  }
  else {							
    t.df <- combn(a, 2, function(ij) {		
      t <- abs(diff(Mean[ij]))/sqrt(sum(Variance[ij]/n[ij]))
      df <- sum(Variance[ij]/n[ij])^2/sum((Variance[ij]/n[ij])^2/(n[ij]-1))
      return(c(t, df))} )
    t <- t.df[1,]
    df <- t.df[2,]
    p <- ptukey(t*sqrt(2), a, df, lower.tail=FALSE)	
    Games.Howell <- cbind(t, df, p)			
    rownames(Games.Howell) <- combn(a, 2, paste, collapse=":")
    return(list(result1=result1, Games.Howell=Games.Howell))
  }
}
# Pairwise Pairwise Multiple Multiple Comparison Comparison Procedures Procedures with Unequal Unequal N‘s
# and/orVariances: A Monte Carlo Study
# Journal Journal of Educational Statistics Educational Statistics,Vol.1, ,Vol.1, No. 2, 1976, pp. 113 . 2, 1976, pp. 113-125

for (i in 3:13){
  jj<-oneway.test(WelAnova[[i]]~Plot,data=WelAnova)
  capture.output(jj,file=paste(Welnames[i],".doc",sep="ANOVA"))
  tt<-tukey(WelAnova[[i]],WelAnova$Plot,method="Games-Howell")Z
  print(tt)
  capture.output(tt,file=paste(Welnames[i],".doc",sep="GHPosthoc"))
}

library(nparcomp)
library(pgirmess) # for KW posthoc comparisons
Gnames<-(Names[1:129,1])
for (i in 80:(ncol(Gradient))){
uu<-kruskalmc(Gradient[[i]]~Plot,data=Gradient) # post-hoc test for the KW test
capture.output(uu,file=paste(Gnames[i],".doc",sep="KW_posthoc"))
}
kruskalmc(PMNM~Plot,data=Gradient) # post-hoc test for the KW test
names(Gradient[80:129])
