#####################################################################################
########################### ENVIRONMENT AND SPECIES #################################
#####################################################################################

# ANALYSIS ON BLACK CHERRY BASAL AREA, WHICH HAS SOME BUT NOT A TON OF ZERO INFLATION
# REMOVING THOSE WITH "UNKNOWN" SOIL TYPES

# LOAD LIBRARIES, SET DIRECTORIES, ETC

  rm(list=ls())
  setwd("~/Dropbox/Duke-Research/casestudies/environment")
  loc <- "~/Dropbox/Duke-Research/casestudies/environment/"

  library(ggplot2)   # for plotting
  library(quantreg)
  library(qrjointEB) # for getting estimates -- my version with Left Censoring
  library(reshape)   # for rearranging dataframes in predict function
  library(coda)      # for ancillary function to coef.qrjoint (overwritten from qrjoint by sourced R code below) 
  library(AER)       # for tobit regression
  library(pROC)      # for AUC and ROC curves
  library(splines)   # for b-splines in Logistic Regression
  library(parallel)  # for running in parallel
  library(dplyr)     # For summarizing easily the simulation results
  library(tidyr)     # For reordering matrix of values
  library(ggplot2)   # For plotting eassily the simulation results
  library(xtable)    # for printing final table
  library(gridExtra)   # for arranging all plots together
  library(coda)      # for convergence checkes, ancillary function to coef.qrjoint (overwritten from qrjoint by sourced R code below) 
  
  theme_update(plot.title = element_text(hjust = 0.5))
  
#####################################################################################
#################################### FUNCTIONS ######################################
#####################################################################################
  
  ########## 
  source("~/Dropbox/Duke-Research/casestudies/environment/0CVDataSplitFunction.R")
  
  ##########  
  # Functions for fitting and updating qrmodels and for evaluating metrics against Portnoy, Logistics
  source("~/Dropbox/Duke-Research/sims/zeroinf/0FittingAndLossFunctions.R")  
  
  ##########  
  # Functions for metric calculations
  source("~/Dropbox/Duke-Research/sims/zeroinf/0ComparisonFunctions.R")  
  
  ########## 
  source("~/Dropbox/Duke-Research/casestudies/environment/0MCMCDiagFunction.R")  
  
  # Function to find tau that corresponds to each (censored) observation, various methods
  modfit.all.qrjoint  <- function(Qarray, y, plot=FALSE, zeroinf=FALSE, xlab="x", ylab="", mcmc="One"){
    # Function to find probability (tau_Y) where each Quantile=Y_i.
    # When zeroinf=TRUE, tau_Y is replaced with tau_Y* drawn from Unif(0, tau_Y) for zero-censored locations.
    # Outputs can be compared to uniform distribution to check model fit.
    # 
    # INPUTS
    # Qarray: An array of predicted quantiles (observation x tau x MCMC iteration) on equally-spaced tau grid
    # y: Original data used to fit models
    # plot: If TRUE, plot to compare to a uniform distribution
    # mcmc: Can be "Summary", "One", "Many" describing how to summarize over posterior draws:
    #       "Summary" takes median of all posterior draws,
    #       "One" grabs one MCMC iteration and uses same for all observations,
    #       "Many" uses a random MCMC iteration for each of n predicted Quantile functions.
    taugrid <- seq(0,1,length=dim(Qarray)[2])
    n <- length(y)
    nmcmc <- dim(Qarray)[3]
    
    if(mcmc=="Summary") MCMC <- apply(Qarray, c(1,2), median)
    if(mcmc=="One")     MCMC <- Qarray[,,sample(1:nmcmc,1)]
    if(mcmc=="Many")    MCMC <- t(apply(Qarray,1,function(f) c(f[,sample(1:nmcmc,1)]) ) ) # random per
    
    tau.Y <- unlist(lapply(1:n,function(f) {approx(x=MCMC[f,], y=taugrid, xout=y[f])$y}))
    if(zeroinf) tau.Y[y==0] <- runif(sum(y==0),0,tau.Y[y==0])
    
    if(plot){
      hist(tau.Y, xlab=xlab, ylab=ylab, prob=T, breaks=max(10, round(n/50)))
      #lines(density(tau.Y, from=0, to=1), main="", xlab=xlab, ylab=ylab)
    }
    invisible(tau.Y)
  }
  
  ###
  # Function to cut a variable into decile bins
  cut_dec <- function(x) cut(x, quantile(x,seq(0,1,length=11)),include.lowest=TRUE)
  
#####################################################################################
####################################### DATA ########################################
#####################################################################################
  
  load(file="fiaClusterDataJim.Rdata")  
  dim(x)          # 1223  12, n X Q design matrix with some covariates and a multilevel factor, soils
  colnames(x)     # "stdage", "temp", "deficit", "topo1", "topo2", "topo3"    The continuous vars 
                  # "entvert", "mol", "others", "spodhist", "ult", "ultkan    The indicator/categorical vars

  # Just for reference... do use dummy coded variables in x
  soil <- as.character(rep("alfinc",times=dim(x)[1]))
  soil[x[,"entvert"]==1] <- "entvert"        # 'Entisols','Vertisols'
  soil[x[,"mol"]==1] <- "mol"                # 'Mollisols’
  soil[x[,"spodhist"]==1] <- "spodhist"      # 'Spodosols','Histosols'
  soil[x[,"ult"]==1] <- "ult"                # 'Ultisols_Udults_Others','Ultisols_Aquults','Ultisols_Others'
  soil[x[,"ultkan"]==1] <- "ultkan"          # 'Ultisols_Udults_Kanhapludults'
  soil[x[,"others"]==1] <- "others"          # Removing these

  # Full list of soil types
  #  c('Alfisols','Inceptisols','Entisols','Histosols','Mollisols',
  #    'Spodosols','Vertisols','Ultisols_Udults_Kanhapludults',
  #    'Ultisols_Udults_Others','Ultisols_Aquults','Ultisols_Others’)
  
  # Picking off a couple of response variables to look at
  responses <- data.frame(nRedMaple    = counts[,"ACRU"],   # 99% of sites have red maple
                          nBlackCherry = counts[,"PRSE2"],  # 83% of sites have black cherry
                          baRedMaple    = ba[,"ACRU"],      # 99% of sites have red maple
                          baBlackCherry = ba[,"PRSE2"])     # 83% of sites have black cherry )
  dat.all <- cbind(responses, x, soil)

  # Subset to remove "others"
  dat <- dat.all[dat.all$soil !="others",]
  dat$soil <- factor(dat$soil)

  dim(dat)  # 1211   17
  n <- nrow(dat)  
  
#####################################################################################
######################################  EDA  ########################################
#####################################################################################
# Exploratory Analysis. See 4ECO R files for fuller exploration
  
  qplot(dat$baBlackCherry, geom="histogram", xlab="Black Cherry Basal Area", ylab="", bins=30)  

  summary(dat$baBlackCherry)
  mean(dat$baBlackCherry==0)
  median(dat$baBlackCherry[dat$baBlackCherry!=0])
    
  ggplot(data=dat, aes(x=stdage, y=baBlackCherry)) + geom_point()  
  ggplot(data=dat, aes(x=temp, y=baBlackCherry)) + geom_point()  
  ggplot(data=dat, aes(x=deficit, y=baBlackCherry)) + geom_point()  
  ggplot(data=dat, aes(x=topo1, y=baBlackCherry)) + geom_point()  
  ggplot(data=dat, aes(x=topo2, y=baBlackCherry)) + geom_point()  
  ggplot(data=dat, aes(x=topo3, y=baBlackCherry)) + geom_point()  
  ggplot(data=dat, aes(x=soil, y=baBlackCherry)) + geom_boxplot()

# How these guys got averaged over tracts is making some difference
# If there is no slope, there will be no aspect (true for 13 clusters)
# But I think with these being averaged, we should probably give it some leeway
# Essentially linear relationship at these gradients
  ggplot(data=dat, aes(x=topo1, y=asin(topo1))) + geom_point()
  ggplot(data=dat, aes(x=topo1)) + geom_histogram()
  ggplot(data=dat, aes(x=topo1, y=baBlackCherry)) + geom_point() 

  ggplot(data=dat, aes(x=topo2)) + geom_histogram()
  ggplot(data=dat, aes(x=topo3)) + geom_histogram()

###########################################################################################
#####################################  MODEL FITTING  #####################################
###########################################################################################  

###########################################################################################
# TAU GRID
  
  tau <- round(seq(0.01, 0.99, by=0.01),2)

###########################################################################################
# LOGISTIC REGRESSION

# Training model
fit.lg <- glm(I(baBlackCherry==0) ~ bs(stdage,5) + bs(temp,5) + bs(deficit,5) + 
                bs(topo1,5) + bs(topo2,5) + bs(topo3,5) + as.factor(soil) ,
              data = dat, family="binomial")

prob0.lg <- predict(fit.lg, type="response")

# LOOKING AT ROC CURVES
plot.roc(dat$baBlackCherry==0, prob0.lg)  

###########################################################################################
# TOBIT REGRESSION

  #######
  # Using AER package: defaults to gaussian and left censored at 0
  fit.tb <- tobit(baBlackCherry ~ stdage + temp + deficit + topo1 + topo2 + topo3 + as.factor(soil), data = dat)
  zapsmall(summary(fit.tb)$coef)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
  # log(scale) shows log of estimate of residual variance, corresponding to value found for scale
  dat$pfit.tb <- pnorm(dat$baBlackCherry, mean=predict(fit.tb), sd=summary(fit.tb)$scale)
  dat$pfit.tb[dat$baBlackCherry==0] <- runif(sum(dat$nBlackCherry==0),0,dat$pfit.tb[dat$baBlackCherry==0]) # better but still not good
  
  
###########################################################################################
# PORTNOY, ZERO-INFLATED (CENSORED) QUANTILE REGRESSION

  response <- survival::Surv(dat$baBlackCherry, (dat$baBlackCherry!=0), type = "left") 
  fit.pt <- crq(response ~ stdage + temp + deficit + topo1 + topo2 + topo3 + as.factor(soil) , data = dat,
                        method="Portnoy", tau=tau)
  Sfit <- summary(fit.pt, tau)  # Does some sort of bootstrap or jackknife (not original)
  betas.pt <- sapply(Sfit, function(lo) return(lo$coefficients[,1]))
  tau.pt <- round(sapply(Sfit, function(lo) return(lo$tau)),2)

  X <- model.matrix(~stdage + temp + deficit + topo1 + topo2 + topo3 + as.factor(soil), data=dat)
  pred.pt<- X %*% betas.pt
  pred.pt.0 <- pmax(pred.pt,0)
  stepf <- apply(pred.pt, 1, function(f){stepfun(tau.pt, c(f[1],f))})
  ra <- rearrange(stepf)

  pred.pt.ra <- do.call("rbind", lapply(ra, function(f){environment(f)$y}))
  pred.pt.ra <- pred.pt.ra[,-c(1,ncol(pred.pt.ra))]
  pred.pt.ra.0 <- pmax(pred.pt.ra,0)
  prob0.pt.ra <- apply(pred.pt.ra, 1, function(f) approx(f, tau.pt, 0, rule=1:2)$y)

  y <- dat$baBlackCherry
  n <- length(y)
  tau.Y <- unlist(lapply(1:n,function(f) {approx(x=pred.pt.ra[f,], y=tau.pt, xout=y[f], rule=1:2)$y}))
  tau.Y[y==0 & !is.na(tau.Y)] <- runif(sum(y==0 & !is.na(tau.Y)),0,tau.Y[y==0& !is.na(tau.Y)])
  mean(is.na(tau.Y))
  hist(tau.Y)
  tau.Y[is.na(tau.Y)] <- runif(sum(is.na(tau.Y)), 0, min(tau.pt))  # not sure if this is fair to do?
  hist(tau.Y)

  # Save "probability fit" for 
  dat$pfit.pt <- tau.Y

###########################################################################################
# JOINT QUANTILE REGRESSION

  # These have fit censored fits in them, coefficients, including the fits below
   #fit.qrj <- qrjoint(baBlackCherry ~ stdage + temp + deficit + topo1 + topo2 + topo3 + as.factor(soil), 
   #                   data=dat, cens=ifelse(dat$baBlackCherry==0,2,0), nsamp = 1000, thin = 20) # burn-in
   #fit.qrj <- update(fit.qrj, nadd=500, append=FALSE) # retain
   #coef.qrj <- coef(fit.qrj, burn.perc = 0, nmc = 500, plot=FALSE)

   #Matrix calculation time per 1e3 iterations = 5.66 
   #Initial lp = -1369.6
   #iter = 2000, lp = -595.411 acpt = 0.14 0.20 0.19 0.18 0.18 0.20 0.19 0.19 0.21 0.18 0.19 0.17 0.11 0.14 0.00 
   #iter = 4000, lp = -563.683 acpt = 0.17 0.14 0.14 0.12 0.16 0.12 0.15 0.15 0.14 0.14 0.15 0.16 0.22 0.15 0.01 
   #iter = 6000, lp = -556.503 acpt = 0.14 0.13 0.12 0.17 0.12 0.17 0.15 0.14 0.14 0.14 0.13 0.15 0.18 0.15 0.17 
   #iter = 8000, lp = -538.211 acpt = 0.16 0.15 0.12 0.15 0.15 0.14 0.13 0.14 0.14 0.15 0.14 0.14 0.17 0.14 0.28 
   #iter = 10000, lp = -526.476 acpt = 0.11 0.14 0.14 0.16 0.14 0.14 0.15 0.14 0.14 0.13 0.14 0.13 0.18 0.12 0.29 
   #iter = 12000, lp = -537.668 acpt = 0.13 0.16 0.17 0.14 0.14 0.15 0.15 0.13 0.14 0.15 0.15 0.13 0.16 0.14 0.25 
   #iter = 14000, lp = -567.343 acpt = 0.15 0.14 0.16 0.14 0.16 0.13 0.15 0.16 0.15 0.15 0.16 0.16 0.16 0.17 0.29 
   #iter = 16000, lp = -560.219 acpt = 0.14 0.14 0.11 0.11 0.14 0.16 0.15 0.15 0.12 0.15 0.13 0.14 0.14 0.14 0.22 
   #iter = 18000, lp = -548.784 acpt = 0.15 0.15 0.14 0.12 0.13 0.15 0.11 0.14 0.15 0.12 0.14 0.14 0.14 0.14 0.18 
   #iter = 20000, lp = -532.278 acpt = 0.13 0.13 0.15 0.17 0.12 0.16 0.16 0.15 0.16 0.14 0.15 0.17 0.16 0.13 0.16 
   #elapsed time: 1284 seconds
   # fit.qrj <- update(fit.qrj, nadd=500, append=FALSE) # retain
   #Initial lp = -532.278
   #iter = 1000, lp = -547.132 acpt = 0.13 0.15 0.17 0.14 0.15 0.13 0.15 0.13 0.17 0.15 0.15 0.14 0.15 0.13 0.15 
   #iter = 2000, lp = -556.386 acpt = 0.14 0.18 0.19 0.16 0.14 0.15 0.17 0.17 0.11 0.17 0.14 0.16 0.14 0.16 0.14 
   #iter = 3000, lp = -544.659 acpt = 0.13 0.17 0.14 0.13 0.14 0.15 0.18 0.13 0.13 0.12 0.14 0.16 0.15 0.10 0.14 
   #iter = 4000, lp = -542.57 acpt = 0.13 0.16 0.14 0.13 0.16 0.13 0.13 0.14 0.15 0.17 0.18 0.13 0.11 0.14 0.11 
   #iter = 5000, lp = -562.787 acpt = 0.13 0.14 0.17 0.18 0.17 0.15 0.16 0.17 0.18 0.13 0.20 0.15 0.13 0.14 0.14 
   #iter = 6000, lp = -519.228 acpt = 0.13 0.10 0.09 0.12 0.15 0.09 0.10 0.11 0.12 0.10 0.13 0.10 0.13 0.10 0.10 
   #iter = 7000, lp = -568.7 acpt = 0.15 0.14 0.16 0.20 0.17 0.14 0.17 0.15 0.14 0.13 0.17 0.14 0.16 0.18 0.16 
   #iter = 8000, lp = -553.262 acpt = 0.16 0.12 0.12 0.13 0.13 0.15 0.13 0.12 0.14 0.11 0.14 0.14 0.14 0.15 0.14 
   #iter = 9000, lp = -535.79 acpt = 0.12 0.13 0.16 0.15 0.12 0.16 0.16 0.13 0.15 0.14 0.12 0.16 0.11 0.14 0.16 
   #iter = 10000, lp = -539.049 acpt = 0.17 0.18 0.14 0.18 0.15 0.17 0.13 0.13 0.11 0.11 0.14 0.14 0.16 0.14 0.17 
   #elapsed time: 652 seconds
   
  #summary(fit.qrj, more = TRUE)   ## WAIC.1 = 1042.77 , WAIC.2 = 1041.68 
  #Qall <- predict(fit.qrj, summarize=FALSE)
  #save(fit.qrj, coef.qrj, Qall , file="QRJ_sub.Rdata")
  

  # These have fit censored fits in them, coefficients, including the fits below
  #fit.nonlin <- qrjoint(baBlackCherry ~ stdage + I(temp*(temp>0)) + I(temp*(temp<=0)) + 
  #deficit + I((topo1-0.05)*(topo1 < 0.05)) + I((topo1-0.05)*(topo1>=0.05))  + topo2 + topo3 + as.factor(soil), 
  #                   data=dat, cens=ifelse(dat$baBlackCherry==0,2,0), nsamp = 1000, thin = 20) # burn-in
  #fit.nonlin <- update(fit.nonlin, nadd=500, append=FALSE) # retain
  #Matrix calculation time per 1e3 iterations = 5.12 
  #Initial lp = -2161.03
  #iter = 2000, lp = -838.82 acpt = 0.14 0.17 0.15 0.14 0.16 0.16 0.14 0.15 0.16 0.16 0.15 0.16 0.15 0.16 0.07 0.14 0.00 
  #iter = 4000, lp = -765.129 acpt = 0.16 0.15 0.17 0.17 0.16 0.14 0.18 0.17 0.17 0.16 0.15 0.16 0.11 0.10 0.13 0.14 0.00 
  #iter = 6000, lp = -700.83 acpt = 0.14 0.13 0.11 0.13 0.11 0.11 0.14 0.13 0.10 0.13 0.13 0.15 0.11 0.12 0.08 0.14 0.02 
  #iter = 8000, lp = -625.989 acpt = 0.12 0.13 0.12 0.13 0.13 0.12 0.17 0.10 0.14 0.17 0.15 0.16 0.09 0.11 0.06 0.11 0.04 
  #iter = 10000, lp = -592.628 acpt = 0.16 0.13 0.16 0.12 0.14 0.13 0.13 0.14 0.10 0.11 0.12 0.12 0.17 0.14 0.16 0.14 0.05 
  #iter = 12000, lp = -559.297 acpt = 0.13 0.11 0.14 0.07 0.13 0.10 0.12 0.16 0.07 0.11 0.10 0.09 0.07 0.13 0.11 0.14 0.05 
  #iter = 14000, lp = -528.077 acpt = 0.17 0.11 0.12 0.14 0.14 0.17 0.12 0.14 0.12 0.16 0.15 0.11 0.16 0.14 0.12 0.13 0.08 
  #iter = 16000, lp = -491.323 acpt = 0.16 0.13 0.11 0.12 0.15 0.14 0.15 0.16 0.19 0.16 0.20 0.17 0.22 0.15 0.22 0.15 0.15 
  #iter = 18000, lp = -492.083 acpt = 0.07 0.12 0.08 0.10 0.09 0.08 0.09 0.11 0.11 0.10 0.11 0.13 0.10 0.12 0.14 0.15 0.13 
  #iter = 20000, lp = -478.491 acpt = 0.13 0.16 0.15 0.19 0.11 0.15 0.11 0.15 0.19 0.11 0.13 0.15 0.22 0.14 0.15 0.14 0.14 
  #elapsed time: 1495 seconds
  #Initial lp = -478.491
  #iter = 1000, lp = -479.415 acpt = 0.15 0.15 0.13 0.14 0.14 0.09 0.15 0.13 0.14 0.11 0.13 0.12 0.17 0.11 0.14 0.14 0.15 
  #iter = 2000, lp = -482.446 acpt = 0.12 0.16 0.19 0.12 0.16 0.11 0.17 0.13 0.19 0.14 0.14 0.14 0.14 0.17 0.15 0.11 0.17 
  #iter = 3000, lp = -472.798 acpt = 0.13 0.11 0.15 0.14 0.10 0.09 0.13 0.11 0.14 0.13 0.12 0.13 0.16 0.12 0.14 0.15 0.14 
  #iter = 4000, lp = -475.834 acpt = 0.24 0.19 0.24 0.21 0.21 0.17 0.23 0.20 0.19 0.23 0.20 0.22 0.24 0.21 0.23 0.21 0.24 
  #iter = 5000, lp = -466.746 acpt = 0.15 0.18 0.15 0.18 0.16 0.13 0.22 0.13 0.17 0.15 0.19 0.20 0.25 0.14 0.14 0.12 0.24 
  #iter = 6000, lp = -483.358 acpt = 0.19 0.15 0.17 0.17 0.17 0.14 0.17 0.18 0.15 0.19 0.19 0.19 0.20 0.16 0.15 0.14 0.24 
  #iter = 7000, lp = -468.836 acpt = 0.20 0.17 0.17 0.20 0.19 0.14 0.14 0.20 0.17 0.17 0.19 0.19 0.18 0.15 0.12 0.17 0.25 
  #iter = 8000, lp = -472.148 acpt = 0.20 0.12 0.17 0.21 0.18 0.18 0.10 0.18 0.13 0.13 0.12 0.17 0.15 0.15 0.11 0.14 0.24 
  #iter = 9000, lp = -462.993 acpt = 0.10 0.08 0.15 0.16 0.13 0.15 0.12 0.14 0.12 0.13 0.08 0.11 0.09 0.09 0.09 0.14 0.23 
  #iter = 10000, lp = -468.561 acpt = 0.12 0.16 0.13 0.18 0.13 0.18 0.11 0.13 0.12 0.13 0.09 0.16 0.15 0.13 0.10 0.15 0.24 
  #elapsed time: 743 seconds
  
  # summary(fit.nonlin, more = TRUE)   ## WAIC.1 = 875.22 , WAIC.2 = 867.69  NOT YET CONVERGED
  
  # OVERWRITING THIS WITH BELOW
  #coef.nonlin <- coef(fit.nonlin, burn.perc = 0, nmc = 500, plot=FALSE)
  #Qnonlin <- predict(fit.nonlin, summarize=FALSE)
  #save(fit.nonlin, coef.nonlin, Qnonlin , file="QRJ_sub_nonlin.Rdata")
  
  
  # Initialize at RQ fits
  fit.nonlin <- qrjoint(baBlackCherry ~ stdage + I(temp*(temp>0)) + I(temp*(temp<=0)) + 
                          deficit + I((topo1-0.05)*(topo1 < 0.05)) + I((topo1-0.05)*(topo1>=0.05))  + topo2 + topo3 + as.factor(soil), 
                        data=dat, cens=ifelse(dat$baBlackCherry==0,2,0), nsamp = 1000, thin = 20, par="RQ") # burn-in
  fit.nonlin <- update(fit.nonlin, nadd=500, append=FALSE)
  fit.nonlin <- update(fit.nonlin, nadd=500, append=FALSE)
  #Matrix calculation time per 1e3 iterations = 12.98 
  #Initial lp = -1654.74
  #iter = 2000, lp = -402.449 acpt = 0.17 0.17 0.17 0.16 0.18 0.17 0.17 0.18 0.18 0.18 0.18 0.18 0.19 0.19 0.15 0.14 0.00 
  #iter = 4000, lp = -368.087 acpt = 0.12 0.16 0.13 0.12 0.09 0.14 0.12 0.15 0.14 0.14 0.13 0.15 0.12 0.13 0.20 0.14 0.01 
  #iter = 6000, lp = -380.608 acpt = 0.12 0.17 0.11 0.13 0.16 0.17 0.16 0.13 0.14 0.15 0.16 0.15 0.14 0.15 0.19 0.14 0.12 
  #iter = 8000, lp = -369.682 acpt = 0.15 0.12 0.13 0.14 0.17 0.16 0.13 0.15 0.16 0.16 0.13 0.14 0.15 0.14 0.16 0.13 0.25 
  #iter = 10000, lp = -363.782 acpt = 0.13 0.13 0.12 0.12 0.15 0.13 0.14 0.14 0.14 0.14 0.13 0.15 0.14 0.13 0.14 0.14 0.28 
  #iter = 12000, lp = -356.095 acpt = 0.14 0.13 0.13 0.13 0.15 0.15 0.12 0.16 0.15 0.16 0.16 0.14 0.13 0.15 0.14 0.14 0.30 
  #iter = 14000, lp = -369.439 acpt = 0.13 0.15 0.13 0.13 0.13 0.13 0.13 0.14 0.14 0.14 0.13 0.15 0.14 0.18 0.15 0.15 0.26 
  #iter = 16000, lp = -385.201 acpt = 0.15 0.15 0.17 0.16 0.16 0.14 0.17 0.17 0.15 0.14 0.16 0.14 0.16 0.13 0.14 0.15 0.29 
  #iter = 18000, lp = -373.551 acpt = 0.13 0.14 0.12 0.13 0.13 0.14 0.13 0.11 0.14 0.11 0.14 0.13 0.14 0.12 0.12 0.13 0.22 
  #iter = 20000, lp = -370.273 acpt = 0.13 0.14 0.16 0.14 0.17 0.17 0.16 0.16 0.15 0.16 0.14 0.14 0.17 0.17 0.15 0.14 0.23 
  #elapsed time: 1525 seconds
  #Initial lp = -370.273
  #iter = 1000, lp = -372.018 acpt = 0.14 0.16 0.14 0.13 0.15 0.16 0.12 0.14 0.15 0.16 0.15 0.16 0.16 0.17 0.15 0.12 0.22 
  #iter = 2000, lp = -357.914 acpt = 0.11 0.11 0.11 0.12 0.11 0.11 0.09 0.11 0.09 0.13 0.13 0.15 0.10 0.12 0.11 0.11 0.16 
  #iter = 3000, lp = -362.251 acpt = 0.16 0.15 0.16 0.15 0.18 0.14 0.16 0.18 0.15 0.19 0.15 0.14 0.16 0.15 0.16 0.16 0.20 
  #iter = 4000, lp = -353.475 acpt = 0.17 0.14 0.12 0.15 0.14 0.16 0.14 0.14 0.14 0.14 0.14 0.14 0.16 0.15 0.16 0.15 0.20 
  #iter = 5000, lp = -386.857 acpt = 0.17 0.17 0.15 0.15 0.16 0.19 0.17 0.18 0.16 0.15 0.15 0.19 0.15 0.16 0.17 0.17 0.21 
  #iter = 6000, lp = -365.601 acpt = 0.17 0.18 0.19 0.17 0.12 0.14 0.17 0.14 0.17 0.17 0.17 0.16 0.17 0.16 0.17 0.14 0.22 
  #iter = 7000, lp = -366.728 acpt = 0.11 0.14 0.14 0.14 0.13 0.14 0.15 0.14 0.13 0.14 0.13 0.13 0.12 0.11 0.14 0.15 0.18 
  #iter = 8000, lp = -407.919 acpt = 0.11 0.17 0.16 0.17 0.17 0.13 0.16 0.13 0.15 0.18 0.15 0.19 0.15 0.17 0.16 0.15 0.19 
  #iter = 9000, lp = -378.112 acpt = 0.13 0.14 0.17 0.18 0.16 0.12 0.14 0.13 0.13 0.13 0.15 0.14 0.13 0.14 0.14 0.12 0.16 
  #iter = 10000, lp = -379.335 acpt = 0.17 0.12 0.12 0.15 0.14 0.12 0.16 0.12 0.14 0.12 0.14 0.13 0.14 0.14 0.14 0.15 0.15 
  #elapsed time: 771 seconds
  #Initial lp = -379.335
  #iter = 1000, lp = -382.275 acpt = 0.17 0.18 0.16 0.16 0.19 0.18 0.19 0.18 0.18 0.17 0.17 0.15 0.17 0.16 0.16 0.15 0.22 
  #iter = 2000, lp = -377.06 acpt = 0.16 0.11 0.08 0.08 0.16 0.12 0.11 0.14 0.14 0.12 0.13 0.16 0.16 0.12 0.15 0.17 0.18 
  #iter = 3000, lp = -379.371 acpt = 0.16 0.12 0.12 0.13 0.12 0.14 0.13 0.12 0.13 0.13 0.15 0.14 0.12 0.14 0.15 0.16 0.18 
  #iter = 4000, lp = -391.386 acpt = 0.18 0.15 0.15 0.14 0.16 0.14 0.12 0.17 0.15 0.17 0.15 0.13 0.16 0.17 0.13 0.14 0.17 
  #iter = 5000, lp = -360.987 acpt = 0.17 0.12 0.13 0.13 0.14 0.13 0.11 0.10 0.12 0.13 0.11 0.13 0.15 0.12 0.13 0.14 0.14 
  #iter = 6000, lp = -394.548 acpt = 0.14 0.17 0.18 0.19 0.18 0.14 0.15 0.16 0.21 0.17 0.15 0.17 0.14 0.15 0.14 0.13 0.15 
  #iter = 7000, lp = -375.748 acpt = 0.14 0.15 0.13 0.11 0.13 0.15 0.15 0.16 0.15 0.14 0.15 0.15 0.14 0.15 0.15 0.16 0.15 
  #iter = 8000, lp = -372.284 acpt = 0.11 0.14 0.14 0.12 0.18 0.14 0.16 0.15 0.15 0.14 0.14 0.11 0.16 0.19 0.12 0.12 0.16 
  #iter = 9000, lp = -375.365 acpt = 0.17 0.17 0.14 0.14 0.15 0.17 0.18 0.15 0.13 0.15 0.16 0.16 0.15 0.14 0.17 0.13 0.15 
  #iter = 10000, lp = -383.797 acpt = 0.16 0.15 0.10 0.12 0.10 0.12 0.11 0.13 0.13 0.16 0.13 0.14 0.14 0.11 0.14 0.13 0.12 
  #elapsed time: 770 seconds
  
  #summary(fit.nonlin, more = TRUE)   # 711.54 , WAIC.2 = 707.38 
  #coef.nonlin <- coef(fit.nonlin, burn.perc = 0, nmc = 500, plot=FALSE)
  #Qnonlin <- predict(fit.nonlin, summarize=FALSE)
  #save(fit.nonlin, coef.nonlin, Qnonlin , file="QRJ_sub_nonlin.Rdata")
  
  # Save probability fits for each model 
  load("QRJ_sub.Rdata")
  dat$pfit.qrj <- modfit.all.qrjoint(Qall, dat$baBlackCherry, plot=F, zeroinf=T, mcmc="Summary") 
  
  load("QRJ_sub_nonlin.Rdata")
  # Save probability fit for 
  dat$pfit.qrj.nonlin <- modfit.all.qrjoint(Qnonlin, dat$baBlackCherry, plot=F, zeroinf=T, mcmc="Summary") 
  dat$pfit.emp <- (rank(dat$baBlackCherry, ties.method="min") - 1)/length(dat$baBlackCherry)
  
  
#####################################################################################
############################# COEFFICIENT PLOTS #####################################
#####################################################################################
# Coefficient plots to compare the qrjoint fit estimates to the Portnoy estimates
    
# LINEAR COVARIATE MODEL
  tau <- round(seq(0.01, 0.99, by=0.01),2)
  load("~/Dropbox/Duke-Research/casestudies/environment/QRJ_sub.Rdata")
  
  fit.gp <- fit.qrj
  rm(fit.qrj,coef.qrj)
  
  beta.hat <- coef(fit.gp, plot = FALSE, burn = 0.1, nmc = 500)   
  beta.gp <- beta.hat$beta.est[2:100,,]
  
  beta.tb <- array(NA, dim(beta.gp), dimnames=dimnames(beta.gp))
  beta.pt <- array(NA, dim(beta.gp), dimnames=dimnames(beta.gp))
  
  # Tobit beta estimates
  sighat <- fit.tb$scale
  med <- qnorm(tau, fit.tb$coef[1], sighat)
  beta.tb[,1,] <- cbind(Lo = med - qnorm(0.975,0,1)*sighat/(sqrt(n)*dnorm(med)),
                        Med= med,
                        Hi = med + qnorm(0.975,0,1)*sighat/(sqrt(n)*dnorm(med)))
  for (i in 2:12){
    beta.tb[,i,] <-  cbind(rep(confint(fit.tb)[i,1], times=length(tau)),
                           rep(fit.tb$coef[i], times=length(tau)),
                           rep(confint(fit.tb)[i,2], times=length(tau)))
  }
  
  # Get fits from Portnoy estimates (place in context of full tau grid)  
  for(j in 1:12) beta.pt[which(tau %in% tau.pt),j,] <- t(sapply(Sfit, function(lo) return(lo[[2]][j,c(2,1,3)]))) # get low, med, high
  
  
  pdf(file = paste(loc,"paper/ExampleMultivarCoefsAll.pdf",sep=""), height = 5.5, width = 7.5)
    par(mfrow = c(3,4), mar = c(2,4,2,1.1))
    varnames <- gsub("as.factor","",dimnames(beta.hat$beta.samp)$beta)
  
  for(i in 1:12){
    # Joint QR
    getBands(beta.hat$beta.samp[,i,], col = 4, xlab="", ylab=bquote(beta ~ .(varnames[i])))
    
    # Portnoy/Censored single quantiles
    lines(tau, beta.pt[,i,2], col=2, lty = 1)
    lines(tau, beta.pt[,i,1], col=2, lty = 1, lwd=.8)
    lines(tau, beta.pt[,i,3], col=2, lty = 1, lwd=.8)
    
    # Tobit Regression
    lines(tau, beta.tb[,i,2], col='gray') 
    lines(tau, beta.tb[,i,1], col='gray', lty=2) 
    lines(tau, beta.tb[,i,3], col='gray', lty=2)
    
    abline(h=0)
    
    if(i==1) {
    legend("topleft", c("Joint","Portnoy","Tobit"), col=c(4,2,'gray'), lty=1,lwd=2)
    }
  }
  dev.off()   
  

# NON-LINEAR COVARIATE MODEL

  # Fits
  load("~/Dropbox/Duke-Research/casestudies/environment/QRJ_sub_nonlin.Rdata")
  fit.gp.nonlin <- fit.nonlin
  rm(fit.nonlin, coef.nonlin)
  
  fit.tb.nonlin <- tobit(baBlackCherry ~ stdage + I(temp*(temp>0)) + I(temp*(temp<=0)) + deficit + I((topo1-0.05)*(topo1 < 0.05)) + I((topo1-0.05)*(topo1>=0.05))  + topo2 + topo3 + as.factor(soil), data = dat)
  response <- survival::Surv(dat$baBlackCherry, (dat$baBlackCherry!=0), type = "left") 
  fit.pt.nonlin <- crq(response ~ stdage + I(temp*(temp>0)) + I(temp*(temp<=0)) + deficit + I((topo1-0.05)*(topo1 < 0.05)) + I((topo1-0.05)*(topo1>=0.05))  + topo2 + topo3 + as.factor(soil) , data = dat,
                       method="Portnoy", tau=tau)
  Sfit.nonlin <- summary(fit.pt.nonlin, tau)  # Does some sort of bootstrap or jackknife (not original)
  tau.pt <- round(sapply(Sfit.nonlin, function(lo) return(lo$tau)),2)
  
  # Coefficient estimates
  beta.hat <- coef(fit.gp.nonlin, plot = FALSE, burn = 0.1, nmc = 500)   
  beta.gp.nonlin <- beta.hat$beta.est[2:100,,]
  beta.tb.nonlin <- array(NA, dim(beta.gp.nonlin), dimnames=dimnames(beta.gp.nonlin))
  beta.pt.nonlin <- array(NA, dim(beta.gp.nonlin), dimnames=dimnames(beta.gp.nonlin))

  sighat <- fit.tb.nonlin$scale
  med <- qnorm(tau, fit.tb.nonlin$coef[1], sighat)
  beta.tb.nonlin[,1,] <- cbind(Lo = med - qnorm(0.975,0,1)*sighat/(sqrt(n)*dnorm(med)),
                      Med= med,
                      Hi = med + qnorm(0.975,0,1)*sighat/(sqrt(n)*dnorm(med)))
  for (i in 2:14){
     beta.tb.nonlin[,i,] <-  cbind(rep(confint(fit.tb.nonlin)[i,1], times=length(tau)),
                         rep(fit.tb.nonlin$coef[i], times=length(tau)),
                         rep(confint(fit.tb.nonlin)[i,2], times=length(tau)))
  }

  # Get fits from Portnoy estimates (place in context of full tau grid)  
  for(j in 1:14) beta.pt.nonlin[which(tau %in% tau.pt),j,] <- t(sapply(Sfit.nonlin, function(lo) return(lo[[2]][j,c(2,1,3)]))) # get low, med, high

  pdf(file = paste(loc,"paper/ExampleMultivarCoefsAllNonlin.pdf",sep=""), height = 7.5, width = 7.5)
  par(mfrow = c(4,4), mar = c(2,4,2,1.1))
  varnames <- gsub("as.factor","",dimnames(beta.hat$beta.samp)$beta)

  for(i in 1:14){
    # Joint QR
    getBands(beta.hat$beta.samp[,i,], col = 4, xlab="", ylab=bquote(beta ~ .(varnames[i])))
  
    # Portnoy/Censored single quantiles
    lines(tau, beta.pt.nonlin[,i,2], col=2, lty = 1)
    lines(tau, beta.pt.nonlin[,i,1], col=2, lty = 1, lwd=.8)
    lines(tau, beta.pt.nonlin[,i,3], col=2, lty = 1, lwd=.8)
  
    # Tobit Regression
    lines(tau, beta.tb.nonlin[,i,2], col='gray') 
    lines(tau, beta.tb.nonlin[,i,1], col='gray', lty=2) 
    lines(tau, beta.tb.nonlin[,i,3], col='gray', lty=2)
  
    abline(h=0)
  
    if(i==1) {
     legend("topleft", c("Joint","Portnoy","Tobit"), col=c(4,2,'gray'), lty=1,lwd=2)
    }
  }
  dev.off()


#####################################################################################
############################# DIAGNOSTIC PLOTS ######################################
#####################################################################################

# Tobit diagnostic plots
ggplot(data=dat, aes(x=pfit.tb)) + geom_histogram(breaks=seq(0,1,length=21)) + xlab("Estimated data response proportion")
dat.qq <- as.data.frame(do.call("cbind",qqplot(dat$pfit.tb, qunif(ppoints(length(dat$pfit.tb))), plot.it=FALSE)))
p1a <- ggplot(data=dat.qq, aes(x=x, y=y)) + geom_point() + xlab("actual") + ylab("theoretical") + ggtitle("Tobit Model") + geom_abline(intercept=0, slope=1)

ggplot(data=dat, aes(x=stdage, y=pfit.tb)) + geom_point() + geom_smooth(se=F) + ylim(0,1) + ylab("Estimated data response proportion")
ggplot(data=dat, aes(x=factor(cut_dec(stdage),labels=1:10), y=pfit.tb)) + geom_violin() + xlab("Decile bins of stdage") + ylab("")

ggplot(data=dat, aes(x=temp, y=pfit.tb)) + geom_point() + geom_smooth(se=F) + ylim(0,1) + ylab("Estimated data response proportion")
ggplot(data=dat, aes(x=factor(cut_dec(temp),labels=1:10), y=pfit.tb)) + geom_violin() + xlab("Decile bins of temp") + ylab("")

ggplot(data=dat, aes(x=deficit, y=pfit.tb)) + geom_point() + geom_smooth(se=F) + ylim(0,1) + ylab("Estimated data response proportion")
ggplot(data=dat, aes(x=factor(cut_dec(deficit),labels=1:10), y=pfit.tb)) + geom_violin() + xlab("Decile bins of deficit") + ylab("")

ggplot(data=dat, aes(x=topo1, y=pfit.tb)) + geom_point() + geom_smooth(se=F) + ylim(0,1) + ylab("Estimated data response proportion")
ggplot(data=dat, aes(x=factor(cut_dec(topo1),labels=1:10), y=pfit.tb)) + geom_violin() + xlab("Decile bins of topo1") + ylab("")

ggplot(data=dat, aes(x=topo2, y=pfit.tb)) + geom_point() + geom_smooth(se=F) + geom_vline(xintercept=0, col="red") + ylim(0,1) + ylab("Estimated data response proportion")
ggplot(data=dat, aes(x=factor(cut_dec(topo2),labels=1:10), y=pfit.tb)) + geom_violin() + xlab("Decile bins of topo2") + ylab("")

ggplot(data=dat, aes(x=topo3, y=pfit.tb)) + geom_point() + geom_smooth(se=F) + geom_vline(xintercept=0, col="red") + ylim(0,1) + ylab("Estimated data response proportion")
ggplot(data=dat, aes(x=factor(cut_dec(topo3),labels=1:10), y=pfit.tb)) + geom_violin() + xlab("Decile bins of topo3") + ylab("")

ggplot(data=dat, aes(x=soil, y=pfit.tb)) + geom_point() + ylab("Estimated data response proportion")
ggplot(data=dat, aes(x=soil, y=pfit.tb)) + geom_violin() + ylab("")


# OVERALL PROBABLITY PLOTS

ggplot(data=dat, aes(x=pfit.qrj)) + geom_histogram(breaks=seq(0,1,length=21)) + xlab("Estimated data response proportion")
dat.qq <- as.data.frame(do.call("cbind",qqplot(dat$pfit.qrj, qunif(ppoints(length(dat$pfit.qrj))), plot.it=FALSE)))
p1b <- ggplot(data=dat.qq, aes(x=x, y=y)) + geom_point() + xlab("actual") + ylab("theoretical") + ggtitle("Joint QR Model") + geom_abline(intercept=0, slope=1)

pdf(paste(loc,"paper/DiagnosticPlotsOverall.pdf",sep=""),width=8,height=3)  
grid.arrange(p1a, p1b, ncol = 2)
dev.off()


# Diagnosing non-linearities

p2a <- ggplot(data=dat, aes(x=stdage, y=pfit.qrj)) + geom_point() + geom_smooth(se=F) + ylim(0,1) + ylab("Estimated data response proportion")
p2b <- ggplot(data=dat, aes(x=factor(cut_dec(stdage),labels=1:10), y=pfit.qrj)) + geom_violin() + xlab("Decile bins of stdage") + ylab("")

p3a <- ggplot(data=dat, aes(x=temp, y=pfit.qrj)) + geom_point() + geom_smooth(se=F) + ylim(0,1) + ylab("Estimated data response proportion")
p3b <- ggplot(data=dat, aes(x=factor(cut_dec(temp),labels=1:10), y=pfit.qrj)) + geom_violin() + xlab("Decile bins of temp") + ylab("")

p4a <- ggplot(data=dat, aes(x=deficit, y=pfit.qrj)) + geom_point() + geom_smooth(se=F) + ylim(0,1) + ylab("Estimated data response proportion")
p4b <- ggplot(data=dat, aes(x=factor(cut_dec(deficit),labels=1:10), y=pfit.qrj)) + geom_violin() + xlab("Decile bins of deficit") + ylab("")

p5a <- ggplot(data=dat, aes(x=topo1, y=pfit.qrj)) + geom_point() + geom_smooth(se=F) + ylim(0,1) + ylab("Estimated data response proportion")
p5b <- ggplot(data=dat, aes(x=factor(cut_dec(topo1),labels=1:10), y=pfit.qrj)) + geom_violin() + xlab("Decile bins of topo1") + ylab("")

p6a <- ggplot(data=dat, aes(x=topo2, y=pfit.qrj)) + geom_point() + geom_smooth(se=F) + geom_vline(xintercept=0, col="red") + ylim(0,1) + ylab("Estimated data response proportion")
p6b <- ggplot(data=dat, aes(x=factor(cut_dec(topo2),labels=1:10), y=pfit.qrj)) + geom_violin() + xlab("Decile bins of topo2") + ylab("")

p7a <- ggplot(data=dat, aes(x=topo3, y=pfit.qrj)) + geom_point() + geom_smooth(se=F) + geom_vline(xintercept=0, col="red") + ylim(0,1) + ylab("Estimated data response proportion")
p7b <- ggplot(data=dat, aes(x=factor(cut_dec(topo3),labels=1:10), y=pfit.qrj)) + geom_violin() + xlab("Decile bins of topo3") + ylab("")

p8a <- ggplot(data=dat, aes(x=soil, y=pfit.qrj)) + geom_point() + ylab("Estimated data response proportion")
p8b <- ggplot(data=dat, aes(x=soil, y=pfit.qrj)) + geom_violin() + ylab("")

pdf(paste(loc,"paper/DiagnosticPlotsQRJointOverall.pdf",sep=""),width=8,height=3)  
grid.arrange(p1a, p1b, ncol = 2)
dev.off()

pdf(paste(loc,"paper/DiagnosticPlotsQRJointLinearity.pdf",sep=""),width=8,height=12)  
grid.arrange(p3a, p3b, p5a, p5b, p4a, p4b, p8a, p8b, ncol = 2)
dev.off()


# Plots from Non-linear or spline model
p1a <- ggplot(data=dat, aes(x=pfit.qrj.nonlin)) + geom_histogram(breaks=seq(0,1,length=21)) + xlab("Estimated data response proportion")
dat.qq <- as.data.frame(do.call("cbind",qqplot(dat$pfit.qrj.nonlin, qunif(ppoints(length(dat$pfit.qrj.nonlin))), plot.it=FALSE)))
p1b <- ggplot(data=dat.qq, aes(x=x, y=y)) + geom_point() + xlab("actual") + ylab("theoretical")

p2a <- ggplot(data=dat, aes(x=stdage, y=pfit.qrj.nonlin)) + geom_point() + geom_smooth(se=F) + ylim(0,1) + ylab("Estimated data response proportion")
p2b <- ggplot(data=dat, aes(x=factor(cut_dec(stdage),labels=1:10), y=pfit.qrj.nonlin)) + geom_violin() + xlab("Decile bins of stdage") + ylab("")

p3a <- ggplot(data=dat, aes(x=temp, y=pfit.qrj.nonlin)) + geom_point() + geom_smooth(se=F) + ylim(0,1) + ylab("Estimated data response proportion")
p3b <- ggplot(data=dat, aes(x=factor(cut_dec(temp),labels=1:10), y=pfit.qrj.nonlin)) + geom_violin() + xlab("Decile bins of temp") + ylab("")

p4a <- ggplot(data=dat, aes(x=deficit, y=pfit.qrj.nonlin)) + geom_point() + geom_smooth(se=F) + ylim(0,1) + ylab("Estimated data response proportion")
p4b <- ggplot(data=dat, aes(x=factor(cut_dec(deficit),labels=1:10), y=pfit.qrj.nonlin)) + geom_violin() + xlab("Decile bins of deficit") + ylab("")

p5a <- ggplot(data=dat, aes(x=topo1, y=pfit.qrj.nonlin)) + geom_point() + geom_smooth(se=F) + ylim(0,1) + ylab("Estimated data response proportion")
p5b <- ggplot(data=dat, aes(x=factor(cut_dec(topo1),labels=1:10), y=pfit.qrj.nonlin)) + geom_violin() + xlab("Decile bins of topo1") + ylab("")

p6a <- ggplot(data=dat, aes(x=topo2, y=pfit.qrj.nonlin)) + geom_point() + geom_smooth(se=F) + geom_vline(xintercept=0, col="red") + ylim(0,1) + ylab("Estimated data response proportion")
p6b <- ggplot(data=dat, aes(x=factor(cut_dec(topo2),labels=1:10), y=pfit.qrj.nonlin)) + geom_violin() + xlab("Decile bins of topo2") + ylab("")

p7a <- ggplot(data=dat, aes(x=topo3, y=pfit.qrj.nonlin)) + geom_point() + geom_smooth(se=F) + geom_vline(xintercept=0, col="red") + ylim(0,1) + ylab("Estimated data response proportion")
p7b <- ggplot(data=dat, aes(x=factor(cut_dec(topo3),labels=1:10), y=pfit.qrj.nonlin)) + geom_violin() + xlab("Decile bins of topo3") + ylab("")

p8a <- ggplot(data=dat, aes(x=soil, y=pfit.qrj.nonlin)) + geom_point() + ylab("Estimated data response proportion")
p8b <- ggplot(data=dat, aes(x=soil, y=pfit.qrj.nonlin)) + geom_violin() + ylab("")


pdf("DiagnosticPlotsQRJointOverallNonlin.pdf",width=8,height=24)  
grid.arrange(p1a, p1b, p2a, p2b, p3a, p3b, p4a, p4b,
             p5a, p5b, p6a, p6b, p7a, p7b, p8a, p8b, ncol = 2)
dev.off()

pdf(paste(loc,"paper/DiagnosticPlotsQRJointSubsetNonlin.pdf",sep=""),width=8,height=6)  
grid.arrange(p3a, p3b, p5a, p5b, ncol = 2)
dev.off()

