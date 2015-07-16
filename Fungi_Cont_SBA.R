#########################
# For SBAPC change. First do all of F_Factor_Plots script
CFun0.1<- glmer(Fungi ~ 1 + (SBAPC| Site), data = Gradient6,family=poisson)
CFun0.2<- glmer(Fungi ~ 1 + (1 | Site), data = GradBD,family=poisson)
AICc(CFun0.1,CFun0.2) # shows that the random effects should include SBAPC change 
CFun1<-glmer(Fungi ~ SBAPC +Dung+(SBAPC| Site), data = Gradient6,family=poisson)
CFun2<-glmer(Fungi ~ SBAPC*Dung + (SBAPC| Site), data = Gradient6,family=poisson)
AICc(CFun0.1,CFun1,CFun2)
anova(CFun0.1,CFun2) # Implies that should keep both. Let's try model selection though
AICc(CFun0.1,CFun1,CFun2) # CFun1 (i.e. the one with the additive term) is better slightly
CFun3<-glmer(Fungi~Dung+ (SBAPC| Site), data = Gradient6,family=poisson)
AICc(CFun0.1,CFun1,CFun2,CFun3)
CFun4<-glmer(Fungi~SBAPC+ (SBAPC| Site), data = Gradient6,family=poisson)
CFun5<-glmer(Fungi~SBAPC+I(SBAPC^2)+ (SBAPC| Site), data = Gradient6,family=poisson)
AICc(CFun0.1,CFun1,CFun2,CFun3,CFun4,CFun5) # So the isolation of the quadratic term is sign and should be left in
# Best model is CFun5, according to AIC. Let's see about the weight of each model
#come up with a list of models 
Modelfun<-list(CFun0.1,CFun1,CFun2,CFun3,CFun4,CFun5)
#summarise these in this table
Model_tab<-model.sel(Modelfun)
Model_tab # Nearly all weight goes to the best model
par(mfrow=c(1,1))
plot(CFun5) # residuals look visually fine
summary(CFun5)
confint(CFun5) 
qqnorm(resid(CFun5))
r.squaredGLMM(CFun0.1);r.squaredGLMM(CFun1); r.squaredGLMM(CFun2);
r.squaredGLMM(CFun3); r.squaredGLMM(CFun4)
r.squaredGLMM(CFun5)

# Use the power of prediction
# random effects of the best model now
Gradient6$Pred_R<-predict(CFun5)
new.data<-expand.grid(SBAPC=seq(0,100,0.1),Site=levels(Gradient6$Site))
new.data$Pred_R<-predict(CFun5,newdata=new.data)
new.data$Pred<-predict(CFun5,newdata=new.data,re.form=NA)
par(mfrow=c(1,1))
theme_set(theme_bw(base_size=12))
Grad_plot1<-ggplot(Gradient6,aes(x=SBAPC,y=Fungi,group=Site,colour=Site))+geom_point()+geom_line(data=new.data,aes(y=exp(Pred_R)))
Grad_plot1
Grad_plot2<-Grad_plot1+geom_line(data=new.data,size=2,colour="black",aes(y=exp(Pred),x=SBA.PC.Diff,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot2+xlab("Percentage loss of basal area relative to reference")+ylab("Fungal richness")

#################################

## Fungi
# Continuous SBA for count data
Mod0.1<- glmer(Fungi ~ 1 + (SBA.PC.Diff| Site), data = Gradient6,family=poisson)
Mod0.2<- glmer(Fungi ~ 1 + (1 | Site), data = Gradient6,family=poisson)
Mod0.3<- glm(Fungi~1,data=Gradient6,family=poisson)
AICc(Mod0.1,Mod0.2,Mod0.3) # shows that the random effects should include SBAPC change 
Mod1<-glmer(Fungi ~ SBAPC+I(SBAPC^2)+ (SBAPC| Site), data = Gradient6,family=poisson)
Mod2<-glmer(Fungi ~ SBAPC+ (SBAPC| Site), data = Gradient6,family=poisson)
Mod3<-glmer(Fungi ~ SBAPC+I(SBAPC^2)+Dung+pH+SBAPC:Dung+ (SBAPC| Site), data = Gradient6,family=poisson)
Mod4<-glmer(Fungi ~ SBAPC+I(SBAPC^2)+Dung+SBAPC:Dung+ (SBAPC| Site), data = Gradient6,family=poisson)
Mod5<-glmer(Fungi ~ SBAPC+I(SBAPC^2)+SBAPC:Dung+ (SBAPC| Site), data = Gradient6,family=poisson)
Mod6<-glmer(Fungi ~ pH+ (SBAPC| Site), data = Gradient6,family=poisson)
Mod7<-glmer(Fungi~Dung+ (SBAPC| Site), data = Gradient6,family=poisson)
Modelfun<-list(Mod1,Mod2,Mod3,Mod4,Mod5,Mod6,Mod7,Mod0.1)
Model_tab<-model.sel(Modelfun)
Model_tab

r.squaredGLMM(Mod1);
r.squaredGLMM(Mod2);
r.squaredGLMM(Mod3);
r.squaredGLMM(Mod4)
r.squaredGLMM(Mod5);
r.squaredGLMM(Mod6)
r.squaredGLMM(Mod7)
r.squaredGLMM(Mod0.1)

par(mfrow=c(2,2))
plot(Mod1) # residuals look visually fine
summary(Mod1) # best model
confint(Mod1)
coefs <- data.frame(coef(summary(Mod1))); coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)));coefs
# Plot graphs from predictions

#calculate confidence intervals for plot
newdat<-expand.grid(SBAPC=seq(0,1,0.01),
                      Site=levels(Gradient6$Site),
                      Dung=mean(Gradient6$Dung),
                      pH=mean(Gradient6$pH),Fungi=0)

mm <- model.matrix(terms(Mod1),newdat)
newdat$Fungi <- predict(Mod1,newdat,re.form=NA)
pvar1 <- diag(mm %*% tcrossprod(vcov(Mod1),mm))
tvar1 <- pvar1+VarCorr(Mod1)$Site[1]
cmult <- 2
newdat <- data.frame(
  newdat
  , plo = newdat$Fungi-cmult*sqrt(pvar1)
  , phi = newdat$Fungi+cmult*sqrt(pvar1)
  , tlo = newdat$Fungi-cmult*sqrt(tvar1)
  , thi = newdat$Fungi+cmult*sqrt(tvar1)
)



theme_set(theme_bw(base_size=12))
Grad_plot1<-ggplot(Gradient6,aes(x=SBAPC*100,y=(Fungi),group=Site,colour=Site))+geom_point()+guides(color = "none")
Grad_plot1
Grad_plot2<-Grad_plot1+geom_line(data=newdat,size=2,colour="black",aes(y=exp(Fungi),x=SBAPC*100,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot3<-Grad_plot2+geom_ribbon(data=newdat,aes(ymax=exp(thi),ymin=exp(tlo)),alpha=0.1,colour=NA)

g2<-Grad_plot3+labs(x="Percentage loss of basal area relative to reference", y="Mycorhizzal Fungi Richness")
g3<-g2+theme(axis.text = element_text(size = 14, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(1), angle = 90),
             axis.title.x = element_text(size = rel(1)))
g4
Fun<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
Fun
ggsave("Fungi_richness.pdf",width = 8,height = 6,units = "in",dpi = 400)



Gradient6$Pred_R<-predict(Mod1)
new.data<-expand.grid(SBAPC=seq(0,1,0.01),
                      Site=levels(Gradient6$Site),
                      Dung=mean(Gradient6$Dung),
                      pH=mean(Gradient6$pH))
new.data$Pred_R<-predict(Mod1,newdata=new.data)
new.data$Pred<-predict(Mod1,newdata=new.data,re.form=NA)


theme_set(theme_bw(base_size=12))
Grad_plot1<-ggplot(Gradient6,aes(x=SBAPC*100,y=(Fungi),group=Site,colour=Site))+geom_point()+guides(color = "none")+geom_line(data=new.data,aes(y=exp(Pred_R)))
Grad_plot1
Grad_plot2<-Grad_plot1+geom_line(data=new.data,size=2,colour="black",aes(y=exp(Pred),x=SBAPC*100,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot2
g2<-Grad_plot2 + theme(axis.text.x=element_text(angle=55, size=30, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=30, vjust=0.5))+
  labs(x="Percentage loss of basal area relative to reference", y="Mycorhizzal Fungi Richness")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(3), angle = 90),
             axis.title.x = element_text(size = rel(2)))
g4
Fun<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
Fun

