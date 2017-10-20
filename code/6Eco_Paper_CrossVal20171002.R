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
############################# CROSS VALIDATION ######################################
#####################################################################################

#################################################
# DATA WITH 1211 OBSERVATIONS, REMOVING "OTHER" CATEGORY

# Make a new data set
set.seed(2001); data.create(nsim=10, subset=(dat$soil !="others"), data.title="paper" )

if(1){  # For a closer look at each training set
  data.title <- "paper"
  for (j in 1:10){
    load(file = paste("kfold/", data.title, formatC(j, width = 2, format = "d", flag = "0"), ".Rdata",sep=""))
    assign(paste("dat.train.",formatC(j, width = 2, format = "d", flag = "0"), sep=""), dat.train)
  }
}

#################################################
# LINEAR DESIGN
set.seed(3002); pp <- mclapply(1:10, do.work, data.title="paper", fit.title= "paper_lin", 
                               design=formula(baBlackCherry ~ stdage + temp + deficit +
                                               topo1  + topo2 + topo3 + as.factor(soil)),
                               nburn=1000, mc.cores = 10)

# Running a second set of fits (start from different location as default is draw from prior)
# on which to do diagnostics
set.seed(3003); pp <- mclapply(1:10, do.work, data.title="paper", fit.title= "paper_lin2", 
                               design=formula(baBlackCherry ~ stdage + temp + deficit +
                                                topo1  + topo2 + topo3 + as.factor(soil)),
                               nburn=2000, mc.cores = 10)

set.seed(3020); pp <- mclapply(1:10, do.loss.work, data.title="paper", fit.title= "paper_lin2", 
                               design=formula(baBlackCherry ~ stdage + temp + deficit + 
                                                topo1 + topo2 + topo3 + as.factor(soil)),
                               mc.cores = 10)
result.plots(nsim=10, fit.title= "paper_lin2")


# Use Diagnostic Function from "OMCMCDiagFunction.R"
pp <- mclapply(1:10, do.diag, fit.title="paper_lin", fit.title2="paper_lin2", nsamp=500)
save(pp, file="diagnostics_lin.Rdata")
load("diagnostics_lin.Rdata")

pp[[1]]   # probably fine
pp[[2]]   # topo1 a little bad
pp[[3]]
pp[[4]]
pp[[5]]  # spod hist and other soil. Bad
pp[[6]]  # fine
pp[[7]]  # soil variables
pp[[8]]  # probably fine, stdage worst
pp[[9]]  # intercept and some soils still doing poorly
pp[[10]] # fine, deficit worst

# Starting at the RQ estimates
set.seed(3004); pp <- mclapply(1:10, do.work, data.title="paper", fit.title= "paper_lin3", 
                               design=formula(baBlackCherry ~ stdage + temp + deficit +
                                                topo1  + topo2 + topo3 + as.factor(soil)), 
                               parstart="RQ",
                               nburn=1000, mc.cores = 10)

set.seed(3020); pp <- mclapply(1:10, do.loss.work, data.title="paper", fit.title= "paper_lin3", 
                               design=formula(baBlackCherry ~ stdage + temp + deficit + 
                                               topo1 + topo2 + topo3 + as.factor(soil)),
                               mc.cores = 10)
result.plots(nsim=10, fit.title= "paper_lin3")



####################
# NON-LINEAR DESIGN WITH topo1 and temp splines
set.seed(2002); pp <- mclapply(1:10, do.work, data.title="paper", fit.title= "paper_linlin", 
                      design=formula(baBlackCherry ~ stdage + I(temp*(temp>0)) + I(temp*(temp<=0)) + deficit +
                                I((topo1-0.05)*(topo1 < 0.05)) + I((topo1-0.05)*(topo1>=0.05))  + topo2 + topo3 + as.factor(soil)),
                      nburn=2500, mc.cores = 10)

# Running a second set of fits (start from different location as default is draw from prior)
# on which to do diagnostics
set.seed(2003); pp <- mclapply(1:10, do.work, data.title="paper", fit.title= "paper_linlin2", 
                      design=formula(baBlackCherry ~ stdage + I(temp*(temp>0)) + I(temp*(temp<=0)) + deficit +
                                I((topo1-0.05)*(topo1 < 0.05)) + I((topo1-0.05)*(topo1>=0.05))  + topo2 + topo3 + as.factor(soil)),
                      nburn=2500, mc.cores = 10)

# Use Diagnostic Function from "OMCMCDiagFunction.R"
pp <- mclapply(1:10, do.diag, fit.title="paper_linlin", fit.title2="paper_linlin2", nsamp=500)
save(pp, file="diagnostics.Rdata")
load("diagnostics.Rdata")

# Updating fits of subsets used for the paper
set.seed(2004); pp <- mclapply(1:10, do.update, fit.title= "paper_linlin", addn=1000, mc.cores = 10)

# Second set of updates
set.seed(2005); pp <- mclapply(1:10, do.update, fit.title= "paper_linlin2", addn=1000, mc.cores = 10)

# Use diagnostics on updated fits
pp <- mclapply(1:10, do.diag, fit.title="paper_linlin", fit.title2="paper_linlin2", nsamp=1000, mc.cores=10)
save(pp, file="diagnostics.Rdata")
load("diagnostics.Rdata")

# Surya thinks of this as an unofficial way to assess convergence. If WAIC 1 is close to WAIC 2 - converged
waic.all <- do.call("rbind",lapply(pp, function(f) return(f[[4]])))
# Close for 1, 3, 4, 5, others are farther off


# Pull off just the beta coefficient scale reduction factors upper tail
psrf.all <- lapply(pp, function(f) return(f[[1]]$psrf[,2]))
psrf.all <- do.call("cbind", psrf.all)
as.matrix(apply(psrf.all,1,median)) 
as.matrix(apply(psrf.all,1,max))  # lower quantiles of upper part of temp, lower and upper topo1 median to upper


# Look at some posterior correlations between the two chains, just for first dataset fit
load(file = paste(loc,"kfold/fit01paper_linlin.Rdata", sep = ""))
o.list <- list(fit.gp)
rm(fit.gp)
load(file = paste("kfold/fit01paper_linlin2.Rdata", sep = ""))
o.list[[2]] <- fit.gp
rm(fit.gp)

beta.samp.list <- mclapply(o.list, function(f) {coef(f, plot = FALSE, burn = 0, nmc = 500)$beta.samp})
beta.list <- mcmc.list(lapply(beta.samp.list, function(bb) return(mcmc(t(rbind(bb[11,,],bb[51,,],bb[91,,]))))))

effectiveSize(beta.list[[1]]) # Some strong autocorrelations topo1 up, & in others

use.these <- c(3,4,17,18,34,35)
cor(beta.list[[1]][,use.these])
plot(as.matrix(beta.list[[1]][,c(2,11)]))

plot(as.matrix(beta.list[[1]][,c(34,35)]))

plot(as.matrix(beta.list[[1]][,c(35)]),as.matrix(beta.list[[2]][,c(35)]))

# posterior means from the two chains
temp <- zapsmall(cbind(apply(as.matrix(beta.list[[1]]),2,mean),apply(as.matrix(beta.list[[2]]),2,mean)))

cbind(temp[1:14,],temp[15:28, ],temp[29:42,])



# Looking at potential scale reduction factors
# Looks good for 1, 3 (beautiful), 4, 5, 8 (could call okay)
pp[[3]]
dat.plot <- dat.train.02
ggplot(data=dat.plot, aes(x=topo1, y=baBlackCherry)) + geom_point()
ggplot(data=dat.plot, aes(x=stdage, y=baBlackCherry)) + geom_point()
ggplot(data=dat.plot, aes(x=temp, y=baBlackCherry)) + geom_point()
pp[[2]] # topo1
pp[[6]] # topo1 a little bad
pp[[7]] # temp and topo1 bad
pp[[9]] # topo1 having trouble
pp[[10]] # temp big problems


# Pull off just the beta coefficient scale reduction factors upper tail
psrf.all <- lapply(pp, function(f) return(f[[1]]$psrf[,2]))
psrf.all <- do.call("cbind", psrf.all)
as.matrix(apply(psrf.all,1,median)) 
as.matrix(apply(psrf.all,1,max))

# RUNNING THIS ALTHOUGH NOT CONVERGED
set.seed(2004); pp <- mclapply(1:10, do.loss.work, data.title="paper", fit.title= "paper_linlin", 
                               design=formula(baBlackCherry ~ stdage + I(temp*(temp>0)) + I(temp*(temp<=0)) + deficit + 
                                                I((topo1-0.05)*(topo1 < 0.05)) + I((topo1-0.05)*(topo1>=0.05))  + topo2 + topo3 + as.factor(soil)),
                               mc.cores = 10)
result.plots(nsim=10, fit.title= "paper_linlin")


# Updating fits of subsets used for the paper
set.seed(2006); pp <- mclapply(c(2,6,7,9,10), do.update, fit.title= "paper_linlin", addn=1000, mc.cores = 5)

# Second set of updates
set.seed(2007); pp <- mclapply(c(2,6,7,9,10), do.update, fit.title= "paper_linlin2", addn=1000, mc.cores = 5)

# What if I only use the sims that I feel like have converged???
result.plots(nsim=10, fit.title= "paper_linlin")


# What if I only use the sims that I feel like have converged???
#result.plots(nsim=NULL, nsub=c(3), fit.title= "paper_linlin")
