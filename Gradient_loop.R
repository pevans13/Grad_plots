#script to loop through different models

rm(list = ls())

#load up packages
library(ggplot2)
library(plyr)
library(lme4)
library(MuMIn)
library(reshape)
install.packages("stargazer")
library(stargazer)

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

## Read functions that will come in useful later
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
  
  # # Plots graphs with standard error
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

#load data
Gradient<- read.csv("FinPlots99.csv")
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


pdf("Figures/Gradient_models.pdf")
#loop to go through different variables and produce models and figures for each of these
str(Gradient)
names(Gradient)
Poivar <- Gradient[c(1:33)] # variables that need a poisson error distribution
Qlovar <- Gradient[c(75:101)] # Variables that are in percentages

for (i in 10:(ncol(Poivar))){
      Modnull1<-glmer(Poivar[[i]]~ 1 +(1|Site),data=Poivar,family="poisson")
      Mod_cont<-glmer(Poivar[[i]] ~ SBAPC + (1| Site), data = Poivar,family="poisson")
      Mod_cont_NL<-glmer(Poivar[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Poivar,family="poisson")
      Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
      print(Model_sel)
      write.csv(Model_sel,paste(names[i],".csv",sep=""),row.names=F)
              }
}

      Preds$Gradient[[5]] <- predict(model.avg(Model_sel),Preds,re.form=NA)
      mm <- model.matrix((Model_sel),Preds)
      pvar1 <- diag(mm %*% tcrossprod(vcov(model.avg(Model_sel))),mm)
      tvar1 <- pvar1+VarCorr(model.avg(Model_sel))$Site[1]

      newdat<-data.frame(SBAPC=seq(0,1,0.01))
      model.matrix(terms(Mod_cont_NL),newdat)
      
      mm <- model.matrix(terms(Mod_cont_NL),newdat)
      newdat$Pred<-predict(Top_model,newdat,re.form=NA)
      colnames(newdat)[2]<-colnames(Gradient[i])
      
      pvar1 <- diag(mm %*% tcrossprod(vcov(Top_model)),mm)
      tvar1 <- pvar1+VarCorr(Top_model)$Site[1]

      cmult <- 2
      newdat <- data.frame(
        Preds
        , plo = newdat$Pred-cmult*sqrt(pvar1)
        , phi = newdat$Pred+cmult*sqrt(pvar1)
        , tlo = newdat$Pred-cmult*sqrt(tvar1)
        , thi = newdat$Pred+cmult*sqrt(tvar1)
      )
      Preds$PredR<-(predict(model.avg(Model_sel),newdata =Preds))
      Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
      theme_set(theme_bw(base_size=12))
      Plot1<-ggplot(data=Gradient,aes(x=SBAPC*100,y=Gradient[[i]],colour=Site,group=Site))+geom_point(size=4,shape=1)
      Plot2<-Plot1+geom_ribbon(data=newdat,aes(ymax=exp(thi),ymin=exp(tlo)),alpha=0.01,colour=NA)
      Plot3<-Plot2+geom_line(data=newdat,size=2,colour="black",aes(y=exp(Gradient[[i]]),x=SBAPC*100,group=NULL))+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))    
      Plot4<-Plot3+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
      Plot5<-Plot4+xlab("Percentage loss of Basal area")+ylab(names[i])+ggtitle(names[i])+ theme(plot.title = element_text(lineheight=.8, face="bold",size=20))
      print(Plot5)
    }
 else{
      Modnull1<-glmer(Gradient[[i]]~ 1 +(1|Site),data=Gradient)
   Mod_cont<-glmer(Gradient[[i]] ~ SBAPC + (1| Site), data = Gradient)
   Mod_cont_NL<-glmer(Gradient[[i]]~SBAPC+I(SBAPC^2)+(1| Site), data = Gradient)
   Model_sel<-model.sel(Modnull1,Mod_cont,Mod_cont_NL,extra = r.squaredGLMM)
   Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradient$Site)
   Preds<-expand.grid(SBAPC=seq(0,1,0.01),Site=Gradient$Site)
   mm <- model.matrix(terms(model.frame(Model_sel),Preds))
   newdat$Fungi <- predict(model.avg(Model_sel),Preds,re.form=NA)
   pvar1 <- diag(mm %*% tcrossprod(vcov(model.avg(Model_sel)),mm))
   tvar1 <- pvar1+VarCorr(model.avg(Model_sel))$Site[1]
   cmult <- 2
   newdat <- data.frame(
     newdat
     , plo = newdat$Gradient[[i]]-cmult*sqrt(pvar1)
     , phi = newdat$Gradient[[i]]+cmult*sqrt(pvar1)
     , tlo = newdat$Gradient[[i]]-cmult*sqrt(tvar1)
     , thi = newdat$Gradient[[i]]+cmult*sqrt(tvar1)
   Preds$PredR<-(predict(model.avg(Model_sel),newdata =Preds))
   Preds$Pred<-(predict(model.avg(Model_sel),newdata =Preds,re.form=NA))
   theme_set(theme_bw(base_size=12))
   Plot1<-ggplot(data=Gradient,aes(x=SBAPC,y=Gradient[[i]],colour=Site,group=Site))+geom_point(size=4,shape=1)
   Plot2<-Plot1+geom_ribbon(data=Preds,aes(ymax=exp(phi),ymin=exp(plo)),alpha=0.01,colour=NA)
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

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
se <- function(x) sqrt(var(x)/length(x))
head(Gradient)

Gradient_melt<-melt(Gradient,id.vars = c("Site","Plot"))

Gradient_summary<-ddply(Gradient_melt,.(Plot,variable),summarise,Mean=mean(value,na.rm = T),SE=se(value))

ggplot(Gradient_summary,aes(x=Plot,y=Mean,ymax=Mean+(2*SE),ymin=Mean-(2*SE)))+geom_pointrange()+geom_line(group=1)+facet_wrap(~variable,scales = "free_y")
ggsave("Figures/Gradient_summary.pdf",height=15,width=30,dpi=500,units="in")
head(Gradient_melt)

theme_set(theme_bw(base_size=12))
Fungi_s<-subset(Gradient_summary, variable=='Fungi')
Plot1<-ggplot(data=Fungi_s,aes(x=Plot,y=Mean))+geom_point()+geom_pointrange(aes(ymin=Mean+SE,ymax=Mean-SE), width = 1)+geom_line(group=1)

