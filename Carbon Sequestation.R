Gradient4<-read.csv("F:/PhD/Results/FinPlots4.csv")
# Get the confidence intervals and se and plot
newSE <- summarySE(Gradient4, measurevar="CarbonSoil", groupvars=c("Plot"))
g<-ggplot(newSE, aes(x=Plot, y=CarbonSoil,group=1)) + 
  geom_errorbar(aes(ymin=CarbonSoil-se, ymax=CarbonSoil+se), width=0.1,size=1.3) +
  geom_line(size=1)+geom_point(size=10,shape=20,col="black")
g
g2<-g + theme(axis.text.x=element_text(angle=55, size=30, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=30, vjust=0.5))+
  labs(x="Stage of collapse", y="Soil carbon (%)")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(3.5), angle = 90, vjust=1.5),
             axis.title.x = element_text(size = rel(3.5)))
g4
bb<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
bb

# CarbonFag
newSE <- summarySE(Gradient4, measurevar="CarbonFag", groupvars=c("Plot"))
g<-ggplot(newSE, aes(x=Plot, y=CarbonFag,group=1)) + 
  geom_errorbar(aes(ymin=CarbonFag-se, ymax=CarbonFag+se), width=0.1,size=1.3) +
  geom_line(size=1)+geom_point(size=10,shape=20,col="black")
g
g2<-g + theme(axis.text.x=element_text(angle=55, size=30, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=30, vjust=0.5))+
  labs(x="Stage of collapse", y="Tree carbon (kg C m-2)")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(3.5), angle = 90, vjust=1.5),
             axis.title.x = element_text(size = rel(3.5)))
g4
cc<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
cc

grid.arrange(bb,cc,ncol=2)

Gradient5<-read.csv("F:/PhD/Results/FinPlots5.csv")
## CarbonFag
hist(Gradient5$CarbonFag)
hist(sqrt(Gradient5$CarbonFag))
hist(log(Gradient5$CarbonFag))
boxplot(Gradient5$CarbonFag)
#Shapiro-Wilk for normaility tests as x has levels (without adjusting for multiple testing). 
do.call("rbind", with(Gradient5, tapply(CarbonFag, Plot,
                                        function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
bartlett.test(resid(lm(CarbonFag~Plot))~Plot,data=Gradient5) # Homogeneity of Variance of residuals
kruskal.test(CarbonFag~Plot,data=Gradient5)
ezANOVA(data=Gradient5, dv=.(tCarbonFag), within=.(CarbonFaglot), wid=.(Site), detailed=TRUE) # Gives info on assumptions and signigifance

## Random effects modelling
Modnull<-lm(tCarbonFag~1,data=Gradient5)
Modnull1<-lmer(tCarbonFag~ 1 +(1|Site),data=Gradient5)
# Test to see if random effects make a difference - judge by std. dev being higher than 0
print(Modnull)
print(Modnull1)
AICc(Modnull,Modnull1)
dotplot(ranef(Modnull1,condVar=TRUE),
        lattice.options=list(layout=c(1,1))) # test to see if intercepts change
# Linear regression
lr1<-lmer(tCarbonFag~CarbonFaglot+(1|Site), data=Gradient5)
summary(lr1)
r.squaredGLMM(lr1)
confint(lr1)
coefs <- data.frame(coef(summary(lr1))); coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value))); coefs
plot(lr1)

# Carbon Soil
Mod0.1<- lmer(CSoil2 ~ 1 + (SBAPC| Site), data = Nutrients)
Mod0.2<- lmer(CSoil2 ~ 1 + (1 | Site), data = Nutrients)
Mod0.3<- lm(CSoil2~1,data=Nutrients)
AICc(Mod0.1,Mod0.2,Mod0.3) 
Mod1<-lmer(CSoil2 ~ SBAPC+I(SBAPC^2) +(SBAPC|Site),  data = Nutrients)
Mod2<-lmer(CSoil2 ~ SBAPC +(SBAPC|Site),  data = Nutrients)
Mod3<-lmer(CSoil2 ~ SBAPC+I(SBAPC^2)+Dung+SBAPC:Dung +(SBAPC|Site),  data = Nutrients)
Mod4<-lmer(CSoil2 ~ SBAPC+I(SBAPC^2)+SBAPC:Dung +(SBAPC|Site),  data = Nutrients)
Mod5<-lmer(CSoil2 ~ pH +(SBAPC|Site),  data = Nutrients)
Mod6<-lmer(CSoil2~Dung +(SBAPC|Site),  data = Nutrients)
Modelfun<-list(Mod1,Mod2,Mod3,Mod4,Mod5,Mod6,Mod0.1)
Model_tab<-model.sel(Modelfun)
Model_tab
r.squaredGLMM(Mod1)

Nutrients$Pred_R<-predict(Mod1)
new.data<-expand.grid(SBAPC=seq(0,1,0.01),
                      Site=levels(Nutrients$Site),
                      Dung=mean(Nutrients$Dung),
                      pH=mean(Nutrients$pH))
new.data$Pred_R<-predict(Mod1,newdata=new.data)
new.data$Pred<-predict(Mod1,newdata=new.data,re.form=NA)
theme_set(theme_bw(base_size=12))
Grad_plot1<-ggplot(Nutrients,aes(x=SBAPC*100,y=plogis(CSoil2),group=Site,colour=Site))+geom_point()+
  guides(color = "none")+geom_line(data=new.data,aes(y=plogis(Pred_R)))
Grad_plot1
Grad_plot2<-Grad_plot1+geom_line(data=new.data,size=2,colour="black",aes(y=plogis(Pred),x=SBAPC*100,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot2
g2<-Grad_plot2 + theme(axis.text.x=element_text(angle=55, size=30, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=30, vjust=0.5))+
  labs(x="Percentage loss of basal area relative to reference", y="Soil carbon (%)")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(3.5), angle = 90),
             axis.title.x = element_text(size = rel(2)))
g4
ccc<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
ccc

# Carbon tree
Mod0.1<- lmer(CarbonFag ~ 1 + (SBAPC| Site), data = Gradient5)
Mod0.2<- lmer(CarbonFag ~ 1 + (1 | Site), data = Gradient5)
Mod0.3<- lm(CarbonFag~1,data=Gradient5)
AICc(Mod0.1,Mod0.2,Mod0.3) 
Mod1<-lmer(CarbonFag ~ SBAPC+I(SBAPC^2) +(SBAPC|Site),  data = Gradient5)
Mod2<-lmer(CarbonFag ~ SBAPC +(SBAPC|Site),  data = Gradient5)
Mod3<-lmer(CarbonFag ~ SBAPC+I(SBAPC^2)+Dung+SBAPC:Dung +(SBAPC|Site),  data = Gradient5)
Mod4<-lmer(CarbonFag ~ SBAPC+I(SBAPC^2)+SBAPC:Dung +(SBAPC|Site),  data = Gradient5)
Mod5<-lmer(CarbonFag ~ pH +(SBAPC|Site),  data = Gradient5)
Mod6<-lmer(CarbonFag~Dung +(SBAPC|Site),  data = Gradient5)
Modelfun<-list(Mod1,Mod2,Mod3,Mod4,Mod5,Mod6,Mod0.1)
Model_tab<-model.sel(Modelfun)
Model_tab
r.squaredGLMM(Mod1)

kruskal.test(CarbonFag~Plot,data=Gradient5) # non-parametric one.way anova equivalent
kruskalmc(CarbonFag ~ Plot, data = Gradient5) # post-hoc test for the KW test

# Linear regression
lr1<-lmer(CarbonFag~Plot+(1|Site), data=Gradient5)
summary(lr1)
r.squaredGLMM(lr1)
confint(lr1)
coefs <- data.frame(coef(summary(lr1))); coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value))); coefs
plot(lr1)
summary(glht(lr1,linfct=mcp(Plot="Tukey")))


Gradient5$Pred_R<-predict(Mod1)
new.data<-expand.grid(SBAPC=seq(0,1,0.01),
                      Site=levels(Gradient5$Site),
                      Dung=mean(Gradient5$Dung),
                      pH=mean(Gradient5$pH))
new.data$Pred_R<-predict(Mod1,newdata=new.data)
new.data$Pred<-predict(Mod1,newdata=new.data,re.form=NA)
theme_set(theme_bw(base_size=12))
Grad_plot1<-ggplot(Gradient5,aes(x=SBAPC*100,y=(CarbonFag),group=Site,colour=Site))+geom_point()+geom_line(data=new.data,aes(y=(Pred_R)))
Grad_plot1
Grad_plot2<-Grad_plot1+geom_line(data=new.data,size=2,colour="black",aes(y=exp(Pred),x=SBAPC*100,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot2
g2<-Grad_plot2 + theme(axis.text.x=element_text(angle=55, size=30, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=30, vjust=0.5))+
  labs(x="Percentage loss of basal area relative to reference", y="Fagus tree carbon (%)")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(3.5), angle = 90),
             axis.title.x = element_text(size = rel(2.5)))
g4
bbb<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA, guides(color = "none")))
bbb

CNRatioSoil
# Get the confidence intervals and se and plot
newSE <- summarySE(Gradient4, measurevar="CNRatioSoil", groupvars=c("Plot"))
g<-ggplot(newSE, aes(x=Plot, y=CNRatioSoil,group=1)) + 
  geom_errorbar(aes(ymin=CNRatioSoil-se, ymax=CNRatioSoil+se), width=0.1,size=1.3) +
  geom_line(size=1)+geom_point(size=10,shape=20,col="black")
g
g2<-g + theme(axis.text.x=element_text(angle=55, size=30, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=30, vjust=0.5))+
  labs(x="Stage of collapse", y="Soil carbon/nitrogen Ratio")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(3), angle = 90),
             axis.title.x = element_text(size = rel(3.5)))
g4
h<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
h
grid.arrange(bbb,ccc,ncol=2)
grid.arrange(h,ncol=2)

# CN Ratio
# Linear regression
lr1<-lmer(CNRatioSoil~Plot+(1|Site), data=Gradient5)
summary(lr1)
r.squaredGLMM(lr1)
confint(lr1)
coefs <- data.frame(coef(summary(lr1))); coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value))); coefs
plot(lr1)
summary(glht(lr1,linfct=mcp(Plot="Tukey")))

Mod0.1<- lmer(CNRatioSoil ~ 1 + (SBAPC| Site), data = Gradient5)
Mod0.2<- lmer(CNRatioSoil ~ 1 + (1 | Site), data = Gradient5)
Mod0.3<- lm(CNRatioSoil~1,data=Gradient5)
AICc(Mod0.1,Mod0.2,Mod0.3) 
Mod1<-lmer(CNRatioSoil ~ SBAPC+I(SBAPC^2) +(1|Site),  data = Gradient5)
Mod2<-lmer(CNRatioSoil ~ SBAPC +(1|Site),  data = Gradient5)
Mod3<-lmer(CNRatioSoil ~ SBAPC+I(SBAPC^2)+Dung+pH+SBAPC:Dung +(1|Site),  data = Gradient5)
Mod4<-lmer(CNRatioSoil ~ SBAPC+I(SBAPC^2)+Dung+SBAPC:Dung +(1|Site),  data = Gradient5)
Mod5<-lmer(CNRatioSoil ~ SBAPC+I(SBAPC^2)+SBAPC:Dung +(1|Site),  data = Gradient5)
Mod6<-lmer(CNRatioSoil ~ pH +(1|Site),  data = Gradient5)
Mod7<-lmer(CNRatioSoil~Dung +(1|Site),  data = Gradient5)
Modelfun<-list(Mod1,Mod2,Mod3,Mod4,Mod5,Mod6,Mod7,Mod0.2)
Model_tab<-model.sel(Modelfun)
Model_tab
r.squaredGLMM(Mod1)

# CN Ratio
Gradient5$Pred_R<-predict(Mod1)
new.data<-expand.grid(SBAPC=seq(0,1,0.01),
                      Site=levels(Gradient5$Site),
                      Dung=mean(Gradient5$Dung),
                      pH=mean(Gradient5$pH))
new.data$Pred_R<-predict(Mod1,newdata=new.data)
new.data$Pred<-predict(Mod1,newdata=new.data,re.form=NA)
theme_set(theme_bw(base_size=12))
Grad_plot1<-ggplot(Gradient5,aes(x=SBAPC*100,y=(CNRatioSoil),group=Site,colour=Site))+geom_point()+guides(color = "none")+geom_line(data=new.data,aes(y=(Pred_R)))
Grad_plot1
Grad_plot2<-Grad_plot1+geom_line(data=new.data,size=2,colour="black",aes(y=(Pred),x=SBAPC*100,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot2
g2<-Grad_plot2 + theme(axis.text.x=element_text(angle=55, size=30, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=30, vjust=0.5))+
  labs(x="Percentage loss of basal area relative to reference", y="Soil C/N ratio")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(3.5), angle = 90),
             axis.title.x = element_text(size = rel(2)))
g4
bc<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
bc

grid.arrange(ccc,bc,ncol=2)

## LOIO
newSE <- summarySE(Gradient6, measurevar="LOIO", groupvars=c("Plot"))
g<-ggplot(newSE, aes(x=Plot, y=LOIO,group=1)) + 
  geom_errorbar(aes(ymin=LOIO-se, ymax=LOIO+se), width=0.1,size=1.3) +
  geom_line(size=1)+geom_point(size=10,shape=20,col="black")
g
g2<-g + theme(axis.text.x=element_text(angle=55, size=30, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=30, vjust=0.5))+
  labs(x="Stage of collapse", y="Loss on ignition (organic)")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(3.5), angle = 90,vjust=1.5),
             axis.title.x = element_text(size = rel(3.5)))
g4
C3<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
C3

## LOIM
newSE <- summarySE(Gradient6, measurevar="LOIM", groupvars=c("Plot"))
g<-ggplot(newSE, aes(x=Plot, y=LOIM,group=1)) + 
  geom_errorbar(aes(ymin=LOIM-se, ymax=LOIM+se), width=0.1,size=1.3) +
  geom_line(size=1)+geom_point(size=10,shape=20,col="black")
g
g2<-g + theme(axis.text.x=element_text(angle=55, size=30, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=30, vjust=0.5))+
  labs(x="Stage of collapse", y="Loss on ignition (mineral)")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(3.5), angle = 90,vjust=1.5),
             axis.title.x = element_text(size = rel(3.5)))
g4
C4<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
C4

grid.arrange(C3,C4,ncol=2)

Gradient6$Bulk_Den
## Bulk Density
newSE <- summarySE(Gradient6, measurevar="Bulk_Den", groupvars=c("Plot"))
g<-ggplot(newSE, aes(x=Plot, y=Bulk_Den,group=1)) + 
  geom_errorbar(aes(ymin=Bulk_Den-se, ymax=Bulk_Den+se), width=0.1,size=1.3) +
  geom_line(size=1)+geom_point(size=10,shape=20,col="black")
g
g2<-g + theme(axis.text.x=element_text(angle=55, size=30, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=30, vjust=0.5))+
  xlab("Stage of collapse")+ ylab(bquote('Bulk density' ~'(g cm^-3*)~')                        
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(3.5), angle = 90, vjust=1.5),
             axis.title.x = element_text(size = rel(3.5)))
g4
bb<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
bb

## Random effects modelling
Modnull<-lm(Bulk_Den~1,data=Gradient6)
Modnull1<-lmer(Bulk_Den~ 1 +(1|Site),data=Gradient6)
Modnull2<- lmer(Bulk_Den~ 1 +(1|Soil_Type),data=Gradient6)
Modnull3<- lmer(Bulk_Den~ 1 +(1|Site)+(1|Soil_Type),data=Gradient6)
# Test to see if random effects make a difference - judge by std. dev being higher than 0
print(Modnull1); print(Modnull2); print (Modnull3); print (Modnull)
AICc(Modnull,Modnull1,Modnull2,Modnull3)
Modelfun<-list(Modnull,Modnull1,Modnull2,Modnull3)
#summarise these in this table
Model_tab<-model.sel(Modelfun)
Model_tab # Nearly all weight goes to the best model

dotplot(ranef(Modnull1,condBulk_Den=TRUE),
        lattice.options=list(layout=c(1,1))) # test to see if intercepts change

# Linear regression
lr1<-lmer(Bulk_Den~Plot+(1|Site), data=Gradient6)
summary(lr1)
r.squaredGLMM(lr1)
confint(lr1)
coefs <- data.frame(coef(summary(lr1)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
plot(lr1)
## Random effects modelling
AICc(Modnull1,Modnull2,Modnull3)

Mod1<- lmer(Bulk_Den~Plot+Dung+(1|Site),data=Gradient6)
Mod2<- lmer(Bulk_Den~Plot+(1|Site),data=Gradient6)
Mod3<- lmer(Bulk_Den~Plot*Dung+(1|Site),data=Gradient6)
Mod4<- lmer(Bulk_Den~Dung+(1|Site),data=Gradient6)
AICc(Mod1, Mod2,Mod3,Mod4,Modnull3)
Modelfun<-list(Mod1,Mod2,Mod3,Mod4,Modnull1)
#summarise these in this table
Model_tab<-model.sel(Modelfun)
Model_tab # Nearly all weight goes to the best model

r.squaredGLMM(Mod1); r.squaredGLMM(Mod2);
r.squaredGLMM(Mod3); r.squaredGLMM(Mod4)
r.squaredGLMM(Modnull3)

# or if glm/lm with no random effects, find out r2 from the code below
R2 <- cor(db$Bulk_Den,predict(Mod1))^2
R2
R2 <- cor(db$Bulk_Den,predict(Mod2))^2
R2
R2 <- cor(db$Bulk_Den,predict(Mod3))^2
R2
R2 <- cor(db$Bulk_Den,predict(Mod4))^2

# Continous SBA for count data
Mod0.1<- lmer(Bulk_Den ~ 1 + (SBAPC| Site), data = Gradient6)
Mod0.2<- lmer(Bulk_Den ~ 1 + (1 | Site), data = Gradient6)
Mod0.3<- lmer(Bulk_Den ~ 1 + (1 | Soil_Type), data = Gradient6 )
Mod0.4<-lm(Bulk_Den~1,data=Gradient6)
AICc(Mod0.1,Mod0.2,Mod0.3,Mod0.4) # shows that the random effects should include SBAPC change 
Modelfun<-list(Mod0.1,Mod0.2,Mod0.3,Mod0.4)
#summarise these in this table
Model_tab<-model.sel(Modelfun)
Model_tab # Nearly all weight goes to the best model

Mod1<-glmer(Bulk_Den ~ SBAPC +Dung+(SBAPC| Site), data = Gradient6,family=poisson)
Mod2<-glmer(Bulk_Den ~ SBAPC*Dung + (SBAPC| Site), data = Gradient6,family=poisson)
AICc(Mod0.1,Mod1,Mod2)
anova(Mod0.1,Mod2) # Implies that should keep both. Let's try model selection though
AICc(Mod0.1,Mod1,Mod2) # Mod1 (i.e. the one with the additive term) is better slightly
Mod3<-glmer(Bulk_Den~Dung+ (SBAPC| Site), data = Gradient6,family=poisson)
AICc(Mod0.1,Mod1,Mod2,Mod3)
Mod4<-glmer(Bulk_Den~SBAPC+ (SBAPC| Site), data = Gradient6,family=poisson)
Mod5<-glmer(Bulk_Den~SBAPC+I(SBAPC^2)+ (SBAPC| Site), data = Gradient6,family=poisson)
Modelfun<-list(Mod1,Mod2,Mod3,Mod4,Mod5,Mod0.1)
Model_tab<-model.sel(Modelfun)
Model_tab # Nearly all weight goes to the best model
r.squaredGLMM(Mod1); r.squaredGLMM(Mod2);r.squaredGLMM(Mod3); r.squaredGLMM(Mod4)
r.squaredGLMM(Mod5); r.squaredGLMM(Mod0.1)
plot(Mod) # residuals look visually fine

# Continous SBA for proportional data
Mod0.1<- lmer(Bulk_Den~ 1 + (SBAPC| Site), data = Gradient6, )
Mod0.2<- lmer(Bulk_Den ~ 1 + (1 | Site), data = Gradient6, )
Mod0.3<- lmer(Bulk_Den ~ 1 + (SBAPC| Soil_Type), data = Gradient6, )
Mod0.4<- lmer(Bulk_Den ~ 1 + (1 | Soil_Type), data = Gradient6, )
Mod0.5<- lmer(Bulk_Den ~ 1 +(1|Site) +(1 | Soil_Type), data = Gradient6, )
print(Mod0.1);print(Mod0.2);print(Mod0.3);print(Mod0.4);print(Mod0.5)
AICc(Mod0.1,Mod0.2,Mod0.3,Mod0.4,Mod0.5) # shows that the random effects should include SBAPC change 
Mod1<-lmer(Bulk_Den ~ SBAPC +Dung+(1|Site), data = Gradient6, )
Mod2<-lmer(Bulk_Den ~ SBAPC*Dung +(1|Site) +(1 | Soil_Type), data = Gradient6, )
Mod3<-lmer(Bulk_Den~Dung+(1|Site) +(1 | Soil_Type), data = Gradient6, )
Mod4<-lmer(Bulk_Den~SBAPC+(1|Site) +(1 | Soil_Type), data = Gradient6, )
Mod5<-lmer(Bulk_Den~SBAPC+I(SBAPC^2)+(1|Site) +(1 | Soil_Type), data = Gradient6,)
AICc(Mod1,Mod2,Mod3,Mod4,Mod5,Mod0.1)
Modelfun<-list(Mod1,Mod2,Mod3,Mod4,Mod5,Mod0.5)
#summarise these in this table
Model_tab<-model.sel(Modelfun)
Model_tab # Nearly all weight goes to the best model
r.squaredGLMM(Mod1); r.squaredGLMM(Mod2);
r.squaredGLMM(Mod3); r.squaredGLMM(Mod4)
r.squaredGLMM(Mod5); r.squaredGLMM(Mod0.1)

summary(Mod1) # best model
confint(Mod1)
# extract coefficients for lme4 package
coefs <- data.frame(coef(summary(Mod)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

# Plot graphs from predictions
Gradient6$Pred_R<-predict(Mod1)
new.data<-expand.grid(SBAPC=seq(0,1,0.01),
                      Dung=mean(Gradient6$Dung),
                      Site=levels(Gradient6$Site))
new.data$Pred_R<-predict(Mod1,newdata=new.data)
new.data$Pred<-predict(Mod1,newdata=new.data,re.form=NA)
theme_set(theme_bw(base_size=12))
Grad_plot1<-ggplot(Gradient6,aes(x=SBAPC*100,y=exp(Bulk_Den),group=Site,colour=Site))+geom_point()+geom_line(data=new.data,aes(y=exp(Pred_R)))
Grad_plot1
Grad_plot2<-Grad_plot1+geom_line(data=new.data,size=2,colour="black",aes(y=exp(Pred),x=SBAPC*100,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot2
Grad_plot2+xlab("Percentage loss of basal area relative to reference")+ylab("Understorey Condition")

