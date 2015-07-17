## Load important packages

library(lme4) # mixed modelling package
library(MuMIn) # For AICc calculation
library(multcomp) # multiple comparisons 
library(SiZer) # Threshold analysis - piecewise regression
install.packages("gridExtra")
library(gridExtra) ## Display graphs in columns and rows (Depends on graphs - sometimes can use ggplots)
library(ggplot2)
library(lattice) # For dotplot (caterpillar plots)
library(pgirmess) # Kruskal-wallis post-hoc

## Read functions that will come in useful later
# Plots graphs with standard error
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

##### Mulitplot - display lots of ggplots at once #####
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


##### Games-Howell ####
# Carries out post-hoc tests for Welch's Anova
tukey <- function(  data,    			# 観察値ベクトル
                    group,					# 群変数ベクトル
                    method=c("Tukey", "Games-Howell"))	# 手法の選択
{
  OK <- complete.cases(data, group)			# 欠損値を持つケースを除く
  data <- data[OK]
  group <- factor(group[OK])
  n <- tapply(data, group, length)			# 各群のケース数
  a <- length(n)						# 群の数
  phi.e <- sum(n)-a					# 誤差分散（群内不偏分散）の自由度
  Mean <- tapply(data, group, mean)			# 各群の平均値
  Variance <- tapply(data, group, var)			# 各群の不偏分散
  result1 <- cbind(n, Mean, Variance)			# 各群の統計量
  rownames(result1) <- paste("Group", 1:a, sep="")
  method <- match.arg(method)
  if (method == "Tukey") {				# Tukey の方法
    v.e <- sum((n-1)*Variance)/phi.e		# 誤差分散（群内不偏分散）
    t <- combn(a, 2, function(ij)			# 対比較
      abs(diff(Mean[ij]))/sqrt(v.e*sum(1/n[ij])) )
    p <- ptukey(t*sqrt(2), a, phi.e, lower.tail=FALSE)	# 有意確率を計算する
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

# Friedman post-hoc test for non-parametric test
friedman.test.with.post.hoc <- function(formu, data, to.print.friedman = T, to.post.hoc.if.signif = T,  to.plot.parallel = T, to.plot.boxplot = T, signif.P = .05, color.blocks.in.cor.plot = T, jitter.Y.in.cor.plot =F)
{
  # formu is a formula of the shape: 	Y ~ X | block
  # data is a long data.frame with three columns:    [[ Y (numeric), X (factor), block (factor) ]]
  
  # Note: This function doesn't handle NA's! In case of NA in Y in one of the blocks, then that entire block should be removed.
  
  
  # Loading needed packages
  if(!require(coin))
  {
    print("You are missing the package 'coin', we will now try to install it...")
    install.packages("coin")
    library(coin)
  }
  
  if(!require(multcomp))
  {
    print("You are missing the package 'multcomp', we will now try to install it...")
    install.packages("multcomp")
    library(multcomp)
  }
  
  if(!require(colorspace))
  {
    print("You are missing the package 'colorspace', we will now try to install it...")
    install.packages("colorspace")
    library(colorspace)
  }
  
  
  # get the names out of the formula
  formu.names <- all.vars(formu)
  Y.name <- formu.names[1]
  X.name <- formu.names[2]
  block.name <- formu.names[3]
  
  if(dim(data)[2] >3) data <- data[,c(Y.name,X.name,block.name)]	# In case we have a "data" data frame with more then the three columns we need. This code will clean it from them...
  
  # Note: the function doesn't handle NA's. In case of NA in one of the block T outcomes, that entire block should be removed.
  
  # stopping in case there is NA in the Y vector
  if(sum(is.na(data[,Y.name])) > 0) stop("Function stopped: This function doesn't handle NA's. In case of NA in Y in one of the blocks, then that entire block should be removed.")
  
  # make sure that the number of factors goes with the actual values present in the data:
  data[,X.name ] <- factor(data[,X.name ])
  data[,block.name ] <- factor(data[,block.name ])
  number.of.X.levels <- length(levels(data[,X.name ]))
  if(number.of.X.levels == 2) { warning(paste("'",X.name,"'", "has only two levels. Consider using paired wilcox.test instead of friedman test"))}
  
  # making the object that will hold the friedman test and the other.
  the.sym.test <- symmetry_test(formu, data = data,	### all pairwise comparisons
                                teststat = "max",
                                xtrafo = function(Y.data) { trafo( Y.data, factor_trafo = function(x) { model.matrix(~ x - 1) %*% t(contrMat(table(x), "Tukey")) } ) },
                                ytrafo = function(Y.data){ trafo(Y.data, numeric_trafo = rank, block = data[,block.name] ) }
  )
  # if(to.print.friedman) { print(the.sym.test) }
  
  
  if(to.post.hoc.if.signif)
  {
    if(pvalue(the.sym.test) < signif.P)
    {
      # the post hoc test
      The.post.hoc.P.values <- pvalue(the.sym.test, method = "single-step")	# this is the post hoc of the friedman test
      
      
      # plotting
      if(to.plot.parallel & to.plot.boxplot)	par(mfrow = c(1,2)) # if we are plotting two plots, let's make sure we'll be able to see both
      
      if(to.plot.parallel)
      {
        X.names <- levels(data[, X.name])
        X.for.plot <- seq_along(X.names)
        plot.xlim <- c(.7 , length(X.for.plot)+.3)	# adding some spacing from both sides of the plot
        
        if(color.blocks.in.cor.plot)
        {
          blocks.col <- rainbow_hcl(length(levels(data[,block.name])))
        } else {
          blocks.col <- 1 # black
        }
        
        data2 <- data
        if(jitter.Y.in.cor.plot) {
          data2[,Y.name] <- jitter(data2[,Y.name])
          par.cor.plot.text <- "Parallel coordinates plot (with Jitter)"
        } else {
          par.cor.plot.text <- "Parallel coordinates plot"
        }
        
        # adding a Parallel coordinates plot
        matplot(as.matrix(reshape(data2,  idvar=X.name, timevar=block.name,
                                  direction="wide")[,-1])  ,
                type = "l",  lty = 1, axes = FALSE, ylab = Y.name,
                xlim = plot.xlim,
                col = blocks.col,
                main = par.cor.plot.text)
        axis(1, at = X.for.plot , labels = X.names) # plot X axis
        axis(2) # plot Y axis
        points(tapply(data[,Y.name], data[,X.name], median) ~ X.for.plot, col = "red",pch = 4, cex = 2, lwd = 5)
      }
      
      if(to.plot.boxplot)
      {
        # first we create a function to create a new Y, by substracting different combinations of X levels from each other.
        subtract.a.from.b <- function(a.b , the.data)
        {
          the.data[,a.b[2]] - the.data[,a.b[1]]
        }
        
        temp.wide <- reshape(data,  idvar=X.name, timevar=block.name,
                             direction="wide") 	#[,-1]
        wide.data <- as.matrix(t(temp.wide[,-1]))
        colnames(wide.data) <- temp.wide[,1]
        
        Y.b.minus.a.combos <- apply(with(data,combn(levels(data[,X.name]), 2)), 2, subtract.a.from.b, the.data =wide.data)
        names.b.minus.a.combos <- apply(with(data,combn(levels(data[,X.name]), 2)), 2, function(a.b) {paste(a.b[2],a.b[1],sep=" - ")})
        
        the.ylim <- range(Y.b.minus.a.combos)
        the.ylim[2] <- the.ylim[2] + max(sd(Y.b.minus.a.combos))	# adding some space for the labels
        is.signif.color <- ifelse(The.post.hoc.P.values < .05 , "green", "grey")
        
        boxplot(Y.b.minus.a.combos,
                names = names.b.minus.a.combos ,
                col = is.signif.color,
                main = "Boxplots (of the differences)",
                ylim = the.ylim
        )
        legend("topright", legend = paste(names.b.minus.a.combos, rep(" ; PostHoc P.value:", number.of.X.levels),round(The.post.hoc.P.values,5)) , fill =  is.signif.color )
        abline(h = 0, col = "blue")
        
      }
      
      list.to.return <- list(Friedman.Test = the.sym.test, PostHoc.Test = The.post.hoc.P.values)
      if(to.print.friedman) {print(list.to.return)}
      return(list.to.return)
      
    }	else {
      print("The results where not significant, There is no need for a post hoc test")
      return(the.sym.test)
    }
  }
}

##
library(lme4)
install.packages("optimx")
require(optimx)   ## for optim optimizers
## (optimx-specific optimizers require explicit gradients --
##  we could use numDeriv::grad, but this seems to defeat
##  the intention)
install.packages("nloptr")
require(nloptr)
install.packages("dfoptim")
require(dfoptim)  ## for nmkb

namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}

## incorporated in lme4 1.1-7
## originally from https://github.com/lme4/lme4/issues/98 :
## nloptWrap <- function(fn, par, lower, upper, control=list(), ...) {
##     defaultControl <- list(xtol_rel = 1e-6, maxeval = 1e5)
##     for (n in names(defaultControl))
##       if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
##     res <- nloptr(x0=par, eval_f=fn, lb=lower, ub=upper, opts=control, ...)
##     ##     ------
##     with(res,list(par=solution,
##                   fval=objective,
##                   feval=iterations,
##                   conv=if (status>0) 0 else status,
##                   message=message))
## }

##' Attempt to re-fit a [g]lmer model with a range of optimizers.
##' The default is to use all known optimizers for R that satisfy the
##' requirements (do not require explicit gradients, allow
##' box constraints), in three categories; (i) built-in
##' (minqa::bobyqa, lme4::Nelder_Mead), (ii) wrapped via optimx
##' (most of optimx's optimizers that allow box constraints require
##' an explicit gradient function to be specified; the two provided
##' here are really base R functions that can be accessed via optimx,
##' (iii) wrapped via nloptr.
##'
##' @param m a fitted model
##' @param meth.tab a matrix (or data.frame) with columns
##' - method  the name of a specific optimization method to pass to the optimizer
##'           (leave blank for built-in optimizers)
##' - optimizer  the \code{optimizer} function to use
##' @param verbose print progress messages?
##' @return a list of fitted \code{merMod} objects
##' @seealso slice, slice2D in the bbmle package
##' @examples
##' library(lme4)
##' gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##'                  data = cbpp, family = binomial)
##' gm_all <- allFit(gm1)
##' t(sapply(gm_all,fixef))              ## extract fixed effects
##' sapply(gm_all,logLik)                ## log-likelihoods
##' sapply(gm_all,getME,"theta")         ## theta parameters
##' !sapply(gm_all,inherits,"try-error") ## was fit OK?
allFit <- function(m, meth.tab = cbind(optimizer=
                                         rep(c("bobyqa","Nelder_Mead", "optimx",  "nloptwrap"),
                                             c(    1,         1,           2,         2)),
                                       method= c("",        "",  "nlminb","L-BFGS-B",
                                                 "NLOPT_LN_NELDERMEAD", "NLOPT_LN_BOBYQA")),
                   verbose=TRUE,
                   
                   maxfun=1e5)
{
  stopifnot(length(dm <- dim(meth.tab)) == 2, dm[1] >= 1, dm[2] >= 2,
            is.character(optimizer <- meth.tab[,"optimizer"]),
            is.character(method    <- meth.tab[,"method"]))
  fit.names <- paste(optimizer, method, sep=".")
  res <- setNames(as.list(fit.names), fit.names)
  for (i in seq_along(fit.names)) {
    if (verbose) cat(fit.names[i],": ")
    ctrl <- list(optimizer=optimizer[i])
    ctrl$optCtrl <- switch(optimizer[i],
                           optimx    = list(method   = method[i]),
                           nloptWrap = list(algorithm= method[i]),
                           list(maxfun=maxfun))
    ctrl <- do.call(if(isGLMM(m)) glmerControl else lmerControl, ctrl)
    tt <- system.time(rr <- tryCatch(update(m, control = ctrl), error = function(e) e))
    attr(rr, "optCtrl") <- ctrl$optCtrl # contains crucial info here
    attr(rr, "time") <- tt  # store timing info
    res[[i]] <- rr
    if (verbose) cat("[OK]\n")
  }
  ## 
  res
}

summary.allfit <- function(object, ...) {
  which.OK <- !sapply(object,is,"error")
  msgs <- lapply(object[which.OK],function(x) x@optinfo$conv$lme4$messages)
  fixef <- t(sapply(object[which.OK],fixef))
  llik <- sapply(object[which.OK],logLik)
  times <- t(sapply(object[which.OK],attr,"time"))
  feval <- sapply(object[which.OK],function(x) x@optinfo$feval)
  sdcor <- t(sapply(object[which.OK],function(x) {
    aa <- as.data.frame(VarCorr(x))
    setNames(aa[,"sdcor"],c(lme4:::tnames(object[which.OK][[1]]),
                            if (isLMM(object[[1]])) "sigma" else NULL))
  }))
  namedList(which.OK,msgs,fixef,llik,sdcor,times,feval)
}
print.summary.allfit <- function(object,...) {
  if (!which.OK==seq(length(object))) {
    cat("some optimizers failed: ",
        paste(names(object)[!which.OK],collapse=","),"\n")
  }
  
}

# deletes NA rows
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


Theme23<-function (base_size = 12, base_family = "") 
{
  theme(
    line = element_line(colour = "black", size = 0.5, linetype = 1, lineend = "butt"), 
    rect = element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1), 
    text = element_text(family = base_family, face = "plain", colour = "black", size = base_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9), 
    
    axis.text = element_text(size = rel(0.8), colour = "grey50"), 
    strip.text = element_text(size = rel(0.8)), 
    axis.line = element_blank(), 
    axis.text.x = element_text(vjust = 1), 
    axis.text.y = element_text(hjust = 1), 
    axis.ticks = element_line(colour = "grey50"), 
    axis.title.x = element_text(), 
    axis.title.y = element_text(angle = 90), 
    axis.ticks.length = unit(0.15, "cm"), 
    axis.ticks.margin = unit(0.1, "cm"), 
    
    legend.background = element_rect(colour = NA), 
    legend.margin = unit(0.2, "cm"), 
    legend.key = element_rect(fill = "grey95", colour = "white"), 
    legend.key.size = unit(1.2, "lines"), 
    legend.key.height = NULL, 
    legend.key.width = NULL, 
    legend.text = element_text(size = rel(0.8)), 
    legend.text.align = NULL, 
    legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0), 
    legend.title.align = NULL, 
    legend.position = "right", 
    legend.direction = NULL, 
    legend.justification = "center", 
    legend.box = NULL, 
    
    panel.background = element_rect(fill = "grey90", colour = NA), 
    panel.border = element_blank(), 
    panel.grid.major = element_line(colour = "white"), 
    panel.grid.minor = element_line(colour = "grey95", size = 0.25), 
    panel.margin = unit(0.25, "lines"), 
    panel.margin.x = NULL, 
    panel.margin.y = NULL, 
    
    strip.background = element_rect(fill = "grey80", colour = NA), 
    strip.text.x = element_text(), 
    strip.text.y = element_text(angle = -90), 
    
    plot.background = element_rect(colour = "white"), 
    plot.title = element_text(size = rel(1.2)), 
    plot.margin = unit(c(1, 1, 0.5, 0.5), "lines"), complete = TRUE)
}

### Read the important file.
Gradient6<-FinPlots6 <- read.csv("F:/PhD/Chapter 1 Gradient Plots/Results/FinPlots6.csv")
View(Gradient6)

### Start the analysis of all the variables for the gradient plots
## The richness of all ground flora (woody and non-woody species) - "CompTot"
# Inspect data
hist(Gradient6$CompTot)
#Shapiro-Wilk for normaility tests as x has levels (without adjusting for multiple testing). 
do.call("rbind", with(Gradient6, tapply(CompTot, Plot,
                                       function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
bartlett.test(resid(lm(CompTot~Plot))~Plot,data=Gradient6) # Homogeneity of Variance of residuals
plot(aov(CompTot~Plot,data=Gradient6)) # diagnostic plots

oneway.test(CompTot~Plot,data=Gradient6) #one-way ANOVA with welch's correction due to heterogeneity of variance
# Following Welch's one.way, Games-Howell post hoc can be used
tukey(Gradient6$CompTot,Gradient6$Plot,method="Games-Howell") # Need a post-hoc test for Welch's correction anova

# Use the summarySE function to determine standard error for data
newSE <- summarySE(Gradient6, measurevar="CompTot", groupvars=c("Plot"))
g<-ggplot(newSE, aes(x=Plot, y=CompTot,group=1)) + 
  geom_errorbar(aes(ymin=CompTot-se, ymax=CompTot+se), width=0.1,size=1.3) +
  geom_line(size=1)+geom_point(size=10,shape=20,col="black")
g
# Change the axis text
g2<-g + theme(axis.text.x=element_text(angle=55, size=14, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=14, vjust=0.5))+
  labs(x="Stage of collapse", y="Total ground flora richness")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(1), angle = 90,vjust=1.5),
             axis.title.x = element_text(size = rel(1)))
g4
# Change the aesthetics
GF1<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
GF1
## Save the figures to file 
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total ground flora richness_SE.pdf",width = 8,height = 6,units = "in",dpi = 400)
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total ground flora richness_SE.jpg",width = 8,height = 6,units = "in",dpi = 400)

ggplot(Gradient6,aes(x=SBAPC,y=CompTot,colour=Site))+geom_point()+facet_wrap(~Site)+geom_smooth(method="glm", family="poisson") # Visualise how the slope differs at each site
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total ground flora richness_Site difference.pdf",width = 8,height = 6,units = "in",dpi = 400)
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total ground flora richness_Site difference.jpg",width = 8,height = 6,units = "in",dpi = 400)

# Linear regression including site as a random effect
lr1<-glmer(CompTot~Plot+(1|Site), data=Gradient6,family=poisson) # Poisson error distribution used because it is count data
summary(lr1)
r.squaredGLMM(lr1)
confint(lr1)
coefs <- data.frame(coef(summary(lr1)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
summary(glht(lr1,linfct=mcp(Plot="Tukey")))

summary(mod<-glmer(CompTot ~ Plot+(1|Site), data = Gradient6,family=poisson))
summary(glht(mod,linfct=mcp(Plot="Tukey")))
 # confint(mod) ## If needed

## Random effects modelling
Modnull<-glm(CompTot~1,data=Gradient6,family=poisson)
Modnull1<-glmer(CompTot~ 1 +(1|Site),data=Gradient6,family=poisson) # Site only as random effect
Modnull2<- glmer(CompTot~ 1 +(1|Site)+(1|Soil_Type),data=Gradient6,family=poisson) # Try site and soil type as random effects
# Test to see if random effects make a difference - judge by std. dev being higher than 0
print(Modnull)
print(Modnull1)
print(Modnull2) # STD for soil is 0 - don't include.
# Use Modnull1

# Plot a dotplot to test to see if intercepts change
Dp<-dotplot(ranef(Modnull1,condVar=TRUE),
        lattice.options=list(layout=c(1,1)))
# Create different models. Focus on plot first just using Plot as a factor
GF1<-glmer(CompTot ~ Plot + (1| Site), data = Gradient6,family=poisson)
# Use mixed models that includes a proxy measure of herbivore pressure. In this case, it is Dung, a measure of dung counts of all herbivores 
Mod1<- glmer(CompTot~Plot+Dung+(1|Site),data=Gradient6,family=poisson)
Mod2<- glmer(CompTot~Plot+(1|Site),data=Gradient6,family=poisson)
Mod3<- glmer(CompTot~Plot*Dung+(1|Site),data=Gradient6,family=poisson)
Mod4<- glmer(CompTot~Dung+(1|Site),data=Gradient6,family=poisson)
AICc(Mod1, Mod2,Mod3,Mod4,Modnull1)
Modelfun<-list(Mod1,Mod2,Mod3,Mod4,Modnull1)
#summarise these in this table
Model_tab<-model.sel(Modelfun)
Model_tab # 0.695 weight for model with plot only - Mod2
# Check diagnostic plots
plot(Mod2) # Looks ok
# Use r-squared below because mixed effects are used
r.squaredGLMM(Mod2)

# Model continous SBA percent change for count data
# Run null models 
Mod0.1<- glmer(CompTot ~1 + (SBAPC| Site), data = Gradient6,family=poisson)
Mod0.2<- glmer(CompTot ~1 + (1 | Site), data = Gradient6,family=poisson)
Mod0.3<-glm(CompTot~1,data=Gradient6,family=poisson)
AICc(Mod0.1,Mod0.2,Mod0.3) # shows that the random effects should include SBAPC change per site
Mod1<-glmer(CompTot ~ SBAPC+Dung + (SBAPC| Site), data = Gradient6,family=poisson)
Mod2<-glmer(CompTot ~ SBAPC*Dung + (SBAPC| Site), data = Gradient6,family=poisson)
Mod3<-glmer(CompTot~Dung+ (SBAPC| Site), data = Gradient6,family=poisson)
Mod4<-glmer(CompTot~SBAPC+ (SBAPC| Site), data = Gradient6,family=poisson)
Mod5<-glmer(CompTot~SBAPC+I(SBAPC^2)+ (SBAPC| Site), data = Gradient6,family=poisson)
## Test which model exhibits most parsimony based on AICc value
AICc(Mod1,Mod2,Mod3,Mod4,Mod5, Mod0.1)
# Best model is Mod5, according to AIC. Let's see about the weight of each model
#come up with a list of models 
ModelGF<-list(Mod1,Mod2,Mod3,Mod4,Mod5, Mod0.1)
#summarise these in this table
Model_tab<-model.sel(ModelGF)
Model_tab # Nearly all weight goes to the best model, which uses a first and second order term only - Mod5
# Check the diagnostic models
plot(Mod5) # Looks fine
r.squaredGLMM(Mod5) # Obtain r2 values
summary(Mod5) # Check significance of terms

# Plot graphs based on the predictions of the best fitting model
Gradient6$Pred_R<-predict(Mod5)
new.data<-expand.grid(SBAPC=seq(0,1,0.01),
                      Dung=mean(Gradient6$Dung),
                      Site=levels(Gradient6$Site))
# Produce new databases from predictions of data
newdat<-expand.grid(SBAPC=seq(0,1,0.01),
                    Site=levels(Gradient6$Site),
                    Dung=mean(Gradient6$Dung),
                    CompTot=0)
# Prodcue confidence intervals
mm <- model.matrix(terms(Mod5),newdat)
newdat$CompTot <- predict(Mod5,newdat,re.form=NA)
pvar1 <- diag(mm %*% tcrossprod(vcov(Mod5),mm))
tvar1 <- pvar1+VarCorr(Mod5)$Site[1]
cmult <- 2
newdat <- data.frame(
  newdat
  , plo = newdat$CompTot-cmult*sqrt(pvar1)
  , phi = newdat$CompTot+cmult*sqrt(pvar1)
  , tlo = newdat$CompTot-cmult*sqrt(tvar1)
  , thi = newdat$CompTot+cmult*sqrt(tvar1)
)
new.data$Pred_R<-predict(Mod5,newdata=new.data)
new.data$Pred<-predict(Mod5,newdata=new.data,re.form=NA) # Doesn't inlcude random effect terms
theme_set(theme_bw(base_size=12))
Grad_plot1<-ggplot(Gradient6,aes(x=SBAPC*100,y=(CompTot),group=Site,colour=Site))+geom_point()+guides(color = "none")
Grad_plot1
Grad_plot2<-Grad_plot1+geom_line(data=newdat,size=2,colour="black",aes(y=exp(CompTot),x=SBAPC*100,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot3<-Grad_plot2+geom_ribbon(data=newdat,aes(ymax=exp(phi),ymin=exp(plo)),alpha=0.01,colour=NA)
g2<-Grad_plot3+labs(x="Percentage loss of basal area relative to reference", y="Total ground flora richness")
g3<-g2+theme(axis.text = element_text(size = 14, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(1), angle = 90),
             axis.title.x = element_text(size = rel(1)))
g4
g5<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
g5
## Save graphs with lower confidence intevals
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total ground flora richness_SBAMod_plo_hi.pdf",width = 8,height = 6,units = "in",dpi = 400)
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total ground flora richness_SBAMod_plo_hi.jpg",width = 8,height = 6,units = "in",dpi = 400)
Grad_plot2<-Grad_plot1+geom_line(data=newdat,size=2,colour="black",aes(y=exp(CompTot),x=SBAPC*100,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot3<-Grad_plot2+geom_ribbon(data=newdat,aes(ymax=exp(thi),ymin=exp(tlo)),alpha=0.01,colour=NA)
g2<-Grad_plot3+labs(x="Percentage loss of basal area relative to reference", y="Total ground flora richness")
g3<-g2+theme(axis.text = element_text(size = 14, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(1), angle = 90),
             axis.title.x = element_text(size = rel(1)))
g4
g5<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
g5
## Save graphs with higher confidence intevals
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total ground flora richness_SBAMod_tlo_hi.pdf",width = 8,height = 6,units = "in",dpi = 400)
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total ground flora richness_SBAMod_tlo_hi.jpg",width = 8,height = 6,units = "in",dpi = 400)

## The richness of all ground flora (woody and non-woody species) - "CompTree"
# Inspect data
hist(Gradient6$CompTree)
#Shapiro-Wilk for normaility tests as x has levels (without adjusting for multiple testing). 
do.call("rbind", with(Gradient6, tapply(CompTree, Plot,
                                        function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
bartlett.test(resid(lm(CompTree~Plot))~Plot,data=Gradient6) # Homogeneity of Variance of residuals
kruskal.test(CompTree~Plot,data=Gradient6) # non-parametric one.way anova equivalent

# Use the summarySE function to determine standard error for data
newSE <- summarySE(Gradient6, measurevar="CompTree", groupvars=c("Plot"))
g<-ggplot(newSE, aes(x=Plot, y=CompTree,group=1)) + 
  geom_errorbar(aes(ymin=CompTree-se, ymax=CompTree+se), width=0.1,size=1.3) +
  geom_line(size=1)+geom_point(size=10,shape=20,col="black")
g
# Change the axis text
g2<-g + theme(axis.text.x=element_text(angle=55, size=14, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=14, vjust=0.5))+
  labs(x="Stage of collapse", y="Total woody species richness")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(1), angle = 90,vjust=1.5),
             axis.title.x = element_text(size = rel(1)))
g4
# Change the aesthetics
GF1<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
GF1
## Save the figures to file 
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total woody species richness_SE.pdf",width = 8,height = 6,units = "in",dpi = 400)
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total woody species richness_SE.jpg",width = 8,height = 6,units = "in",dpi = 400)

ggplot(Gradient6,aes(x=SBAPC,y=CompTree,colour=Site))+geom_point()+facet_wrap(~Site)+geom_smooth(method="glm", family="poisson") # Visualise how the slope differs at each site
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total woody species richness_Site difference.pdf",width = 8,height = 6,units = "in",dpi = 400)
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total woody species richness_Site difference.jpg",width = 8,height = 6,units = "in",dpi = 400)

# Linear regression including site as a random effect
lr1<-glmer(CompTree~Plot+(1|Site), data=Gradient6,family=poisson) # Poisson error distribution used because it is count data
summary(lr1)

## Random effects modelling
Modnull<-glm(CompTree~1,data=Gradient6,family=poisson)
Modnull1<-glmer(CompTree~ 1 +(1|Site),data=Gradient6,family=poisson) # Site only as random effect
Modnull2<- glmer(CompTree~ 1 +(1|Site)+(1|Soil_Type),data=Gradient6,family=poisson) # Try site and soil type as random effects
# Test to see if random effects make a difference - judge by std. dev being higher than 0
print(Modnull)
print(Modnull1)# STD for site is 0 - don't include.
# Use Modnull

# Create different models. Focus on plot first just using Plot as a factor
Mod<-glm(CompTree ~ Plot, data = Gradient6,family=poisson)
# Use mixed models that includes a proxy measure of herbivore pressure. In this case, it is Dung, a measure of dung counts of all herbivores 
Mod1<- glm(CompTree~Plot+Dung ,data=Gradient6,family=poisson)
Mod2<- glm(CompTree~Plot ,data=Gradient6,family=poisson)
Mod3<- glm(CompTree~Plot*Dung ,data=Gradient6,family=poisson)
Mod4<- glm(CompTree~Dung ,data=Gradient6,family=poisson)
AICc(Mod1, Mod2,Mod3,Mod4,Modnull)
Modelfun<-list(Mod1,Mod2,Mod3,Mod4,Modnull)
#summarise these in this table
Model_tab<-model.sel(Modelfun)
Model_tab # 0.695 weight for model with plot only - Mod2
# Check diagnostic plots
plot(Modnull) # Looks ok
# Use r-squared below because mixed effects are used
r.squaredGLMM(Modnull)

# Model continous SBA percent change for count data
# Run null models 
Mod0.1<- glmer(CompTree ~1 + (SBAPC| Site), data = Gradient6,family=poisson)
Mod0.2<- glmer(CompTree ~1 + (1 | Site), data = Gradient6,family=poisson)
Mod0.3<-glm(CompTree~1,data=Gradient6,family=poisson)
AICc(Mod0.1,Mod0.2,Mod0.3) # shows that the random effects should include SBAPC change per site
Mod1<-glm(CompTree ~ SBAPC+Dung  , data = Gradient6,family=poisson)
Mod2<-glm(CompTree ~ SBAPC*Dung  , data = Gradient6,family=poisson)
Mod3<-glm(CompTree~Dung , data = Gradient6,family=poisson)
Mod4<-glm(CompTree~SBAPC , data = Gradient6,family=poisson)
Mod5<-glm(CompTree~SBAPC+I(SBAPC^2) , data = Gradient6,family=poisson)
## Test which model exhibits most parsimony based on AICc value
AICc(Mod1,Mod2,Mod3,Mod4,Mod5, Mod0.3)
# Best model is Mod4, according to AIC. Let's see about the weight of each model
#come up with a list of models 
ModelGF<-list(Mod1,Mod2,Mod3,Mod4,Mod5, Mod0.3)
#summarise these in this table
Model_tab<-model.sel(ModelGF)
Model_tab # Nearly all weight goes to the best model, which uses a first and second order term only - Mod5
# Check the diagnostic models
plot(Mod4) # Looks fine
r.squaredGLMM(Mod4) # Obtain r2 values
summary(Mod4) # Check significance of terms

# Plot graphs based on the predictions of the best fitting model
Gradient6$Pred_R<-predict(Mod4)
new.data<-expand.grid(SBAPC=seq(0,1,0.01),
                      Dung=mean(Gradient6$Dung),
                      Site=levels(Gradient6$Site))
# Produce new databases from predictions of data
newdat<-expand.grid(SBAPC=seq(0,1,0.01),
                    Site=levels(Gradient6$Site),
                    Dung=mean(Gradient6$Dung),
                    CompTree=0)
# Prodcue confidence intervals
mm <- model.matrix(terms(Mod4),newdat)
newdat$CompTree <- predict(Mod4,newdat,re.form=NA)
pvar1 <- diag(mm %*% tcrossprod(vcov(Mod4),mm))
tvar1 <- pvar1+VarCorr(Mod4)$Site[1]
cmult <- 2
newdat <- data.frame(
  newdat
  , plo = newdat$CompTree-cmult*sqrt(pvar1)
  , phi = newdat$CompTree+cmult*sqrt(pvar1)
  , tlo = newdat$CompTree-cmult*sqrt(tvar1)
  , thi = newdat$CompTree+cmult*sqrt(tvar1)
)
new.data$Pred_R<-predict(Mod4,newdata=new.data)
new.data$Pred<-predict(Mod4,newdata=new.data,re.form=NA) # Doesn't inlcude random effect terms
theme_set(theme_bw(base_size=12))
Grad_plot1<-ggplot(Gradient6,aes(x=SBAPC*100,y=(CompTree),group=Site,colour=Site))+geom_point()+guides(color = "none")
Grad_plot1
Grad_plot2<-Grad_plot1+geom_line(data=newdat,size=2,colour="black",aes(y=exp(CompTree),x=SBAPC*100,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot3<-Grad_plot2+geom_ribbon(data=newdat,aes(ymax=exp(phi),ymin=exp(plo)),alpha=0.01,colour=NA)
g2<-Grad_plot3+labs(x="Percentage loss of basal area relative to reference", y="Total woody species richness")
g3<-g2+theme(axis.text = element_text(size = 14, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(1), angle = 90),
             axis.title.x = element_text(size = rel(1)))
g4
g5<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
g5
## Save graphs with lower confidence intevals
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total woody species richness_SBAMod_plo_hi.pdf",width = 8,height = 6,units = "in",dpi = 400)
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total woody species richness_SBAMod_plo_hi.jpg",width = 8,height = 6,units = "in",dpi = 400)
Grad_plot2<-Grad_plot1+geom_line(data=newdat,size=2,colour="black",aes(y=exp(CompTree),x=SBAPC*100,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot3<-Grad_plot2+geom_ribbon(data=newdat,aes(ymax=exp(thi),ymin=exp(tlo)),alpha=0.01,colour=NA)
g2<-Grad_plot3+labs(x="Percentage loss of basal area relative to reference", y="Total woody species richness")
g3<-g2+theme(axis.text = element_text(size = 14, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(1), angle = 90),
             axis.title.x = element_text(size = rel(1)))
g4
g5<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
g5
## Save graphs with higher confidence intevals
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total woody species richness_SBAMod_tlo_hi.pdf",width = 8,height = 6,units = "in",dpi = 400)
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total woody species richness_SBAMod_tlo_hi.jpg",width = 8,height = 6,units = "in",dpi = 400)

# Non-woody species richness
## The richness of all ground flora (woody and non-woody species) - "CompGF"
# Inspect data
hist(Gradient6$CompGF)
#Shapiro-Wilk for normaility tests as x has levels (without adjusting for multiple testing). 
do.call("rbind", with(Gradient6, tapply(CompGF, Plot,
                                        function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) 
bartlett.test(resid(lm(CompGF~Plot))~Plot,data=Gradient6) # Homogeneity of Variance of residuals
oneway.test(CompGF~Plot,data=Gradient6) #one-way ANOVA with welch's correction due to heterogeneity of variance
# Following Welch's one.way, Games-Howell post hoc can be used
tukey(Gradient6$CompGF,Gradient6$Plot,method="Games-Howell") # Need a post-hoc test for Welch's correction anova

# Use the summarySE function to determine standard error for data
newSE <- summarySE(Gradient6, measurevar="CompGF", groupvars=c("Plot"))
g<-ggplot(newSE, aes(x=Plot, y=CompGF,group=1)) + 
  geom_errorbar(aes(ymin=CompGF-se, ymax=CompGF+se), width=0.1,size=1.3) +
  geom_line(size=1)+geom_point(size=10,shape=20,col="black")
g
# Change the axis text
g2<-g + theme(axis.text.x=element_text(angle=55, size=14, vjust=0.5)) + theme(axis.text.y=element_text(angle=0, size=14, vjust=0.5))+
  labs(x="Stage of collapse", y="Total non-woody species richness")
g3<-g2+theme(axis.text = element_text(size = 50, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(1), angle = 90,vjust=1.5),
             axis.title.x = element_text(size = rel(1)))
g4
# Change the aesthetics
GF1<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
GF1
## Save the figures to file 
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total non-woody species richness_SE.pdf",width = 8,height = 6,units = "in",dpi = 400)
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total non-woody species richness_SE.jpg",width = 8,height = 6,units = "in",dpi = 400)

ggplot(Gradient6,aes(x=SBAPC,y=CompGF,colour=Site))+geom_point()+facet_wrap(~Site)+geom_smooth(method="glm", family="poisson") # Visualise how the slope differs at each site
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total non-woody species richness_Site difference.pdf",width = 8,height = 6,units = "in",dpi = 400)
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total non-woody species richness_Site difference.jpg",width = 8,height = 6,units = "in",dpi = 400)

# Linear regression including site as a random effect
lr1<-glmer(CompGF~Plot+(1|Site), data=Gradient6,family=poisson) # Poisson error distribution used because it is count data
summary(lr1)
r.squaredGLMM(lr1)
confint(lr1)
coefs <- data.frame(coef(summary(lr1)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
summary(glht(lr1,linfct=mcp(Plot="Tukey")))

## Random effects modelling
Modnull<-glm(CompGF~1,data=Gradient6,family=poisson)
Modnull1<-glmer(CompGF~ 1 +(1|Site),data=Gradient6,family=poisson) # Site only as random effect
Modnull2<- glmer(CompGF~ 1 +(1|Site)+(1|Soil_Type),data=Gradient6,family=poisson) # Try site and soil type as random effects
# Test to see if random effects make a difference - judge by std. dev being higher than 0
print(Modnull)
print(Modnull1)
print(Modnull2)# STD for Soil is 0 - don't include.
# Use Modnull1

# Create different models. Focus on plot first just using Plot as a factor
Mod<-glmer(CompGF ~ Plot+(1|Site), data = Gradient6,family=poisson)
# Use mixed models that includes a proxy measure of herbivore pressure. In this case, it is Dung, a measure of dung counts of all herbivores 
Mod1<- glmer(CompGF~Plot+Dung +(1|Site),data=Gradient6,family=poisson)
Mod2<- glmer(CompGF~Plot +(1|Site),data=Gradient6,family=poisson)
Mod3<- glmer(CompGF~Plot*Dung +(1|Site),data=Gradient6,family=poisson)
Mod4<- glmer(CompGF~Dung +(1|Site),data=Gradient6,family=poisson)
AICc(Mod1, Mod2,Mod3,Mod4,Modnull1)
Modelfun<-list(Mod1,Mod2,Mod3,Mod4,Modnull1)
#summarise these in this table
Model_tab<-model.sel(Modelfun)
Model_tab # 0.714 weight for model with plot only - Mod2
# Check diagnostic plots
plot(Mod2) # Looks ok
# Use r-squared below because mixed effects are used
r.squaredGLMM(Mod2)

# Model continous SBA percent change for count data
# Run null models 
Mod0.1<- glmer(CompGF ~1 + (SBAPC| Site), data = Gradient6,family=poisson)
Mod0.2<- glmer(CompGF ~1 + (1 | Site), data = Gradient6,family=poisson)
Mod0.3<-glm(CompGF~1,data=Gradient6,family=poisson)
AICc(Mod0.1,Mod0.2,Mod0.3) # shows that the random effects should include SBAPC change per site
Mod1<-glmer(CompGF ~ SBAPC+Dung  + (SBAPC| Site), data = Gradient6,family=poisson)
Mod2<-glmer(CompGF ~ SBAPC*Dung  + (SBAPC| Site), data = Gradient6,family=poisson)
Mod3<-glmer(CompGF~Dung + (SBAPC| Site), data = Gradient6,family=poisson)
Mod4<-glmer(CompGF~SBAPC+ (SBAPC| Site) , data = Gradient6,family=poisson)
Mod5<-glmer(CompGF~SBAPC+I(SBAPC^2)+ (SBAPC| Site) , data = Gradient6,family=poisson)
## Test which model exhibits most parsimony based on AICc value
AICc(Mod1,Mod2,Mod3,Mod4,Mod5, Mod0.1)
#come up with a list of models 
ModelGF<-list(Mod1,Mod2,Mod3,Mod4,Mod5, Mod0.1)
#summarise these in this table
Model_tab<-model.sel(ModelGF)
Model_tab # Nearly all weight goes to the best model, which uses a first and second order term only - Mod5
# Check the diagnostic models
plot(Mod5) # Looks fine
r.squaredGLMM(Mod5) # Obtain r2 values
summary(Mod4) # Check significance of terms

# Plot graphs based on the predictions of the best fitting model
Gradient6$Pred_R<-predict(Mod4)
new.data<-expand.grid(SBAPC=seq(0,1,0.01),
                      Dung=mean(Gradient6$Dung),
                      Site=levels(Gradient6$Site))
# Produce new databases from predictions of data
newdat<-expand.grid(SBAPC=seq(0,1,0.01),
                    Site=levels(Gradient6$Site),
                    Dung=mean(Gradient6$Dung),
                    CompGF=0)
# Prodcue confidence intervals
mm <- model.matrix(terms(Mod4),newdat)
newdat$CompGF <- predict(Mod4,newdat,re.form=NA)
pvar1 <- diag(mm %*% tcrossprod(vcov(Mod4),mm))
tvar1 <- pvar1+VarCorr(Mod4)$Site[1]
cmult <- 2
newdat <- data.frame(
  newdat
  , plo = newdat$CompGF-cmult*sqrt(pvar1)
  , phi = newdat$CompGF+cmult*sqrt(pvar1)
  , tlo = newdat$CompGF-cmult*sqrt(tvar1)
  , thi = newdat$CompGF+cmult*sqrt(tvar1)
)
new.data$Pred_R<-predict(Mod4,newdata=new.data)
new.data$Pred<-predict(Mod4,newdata=new.data,re.form=NA) # Doesn't inlcude random effect terms
theme_set(theme_bw(base_size=12))
Grad_plot1<-ggplot(Gradient6,aes(x=SBAPC*100,y=(CompGF),group=Site,colour=Site))+geom_point()+guides(color = "none")
Grad_plot1
Grad_plot2<-Grad_plot1+geom_line(data=newdat,size=2,colour="black",aes(y=exp(CompGF),x=SBAPC*100,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot3<-Grad_plot2+geom_ribbon(data=newdat,aes(ymax=exp(phi),ymin=exp(plo)),alpha=0.01,colour=NA)
g2<-Grad_plot3+labs(x="Percentage loss of basal area relative to reference", y="Total non-woody species richness")
g3<-g2+theme(axis.text = element_text(size = 14, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(1), angle = 90),
             axis.title.x = element_text(size = rel(1)))
g4
g5<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
g5
## Save graphs with lower confidence intevals
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total non-woody species richness_SBAMod_plo_hi.pdf",width = 8,height = 6,units = "in",dpi = 400)
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total non-woody species richness_SBAMod_plo_hi.jpg",width = 8,height = 6,units = "in",dpi = 400)
Grad_plot2<-Grad_plot1+geom_line(data=newdat,size=2,colour="black",aes(y=exp(CompGF),x=SBAPC*100,group=NULL))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Grad_plot3<-Grad_plot2+geom_ribbon(data=newdat,aes(ymax=exp(thi),ymin=exp(tlo)),alpha=0.01,colour=NA)
g2<-Grad_plot3+labs(x="Percentage loss of basal area relative to reference", y="Total non-woody species richness")
g3<-g2+theme(axis.text = element_text(size = 14, colour = "black"), panel.background = element_rect(fill = "white", colour = NA))
g4<-g3+theme(axis.title.y = element_text(size = rel(1), angle = 90),
             axis.title.x = element_text(size = rel(1)))
g4
g5<-g4+theme(panel.border = element_rect(color="darkred", size=0.5, linetype="solid",fill=NA))
g5
## Save graphs with higher confidence intevals
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total non-woody species richness_SBAMod_tlo_hi.pdf",width = 8,height = 6,units = "in",dpi = 400)
ggsave("F:/PhD/Chapter 1 Gradient Plots/Figures/Total non-woody species richness_SBAMod_tlo_hi.jpg",width = 8,height = 6,units = "in",dpi = 400)
# Delete this