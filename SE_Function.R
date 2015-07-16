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
tukey <- function(  data,					# 観察値ベクトル
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