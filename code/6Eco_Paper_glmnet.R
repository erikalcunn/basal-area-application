#####################################################################################
########################### ENVIRONMENT AND SPECIES #################################
#####################################################################################
# LOAD DATA AND LIBRARIES

rm(list=ls())
setwd("~/Dropbox/Duke-Research/casestudies/environment")
loc <- "~/Dropbox/Duke-Research/casestudies/environment/"

library(glmnet)
library(splines)
library(pROC)      # for AUC and ROC curves
library(parallel)

#####################################################################################
#################################     DATA       ####################################
#####################################################################################

load(file="fiaClusterDataJim.Rdata")
soil <- as.character(rep("alfinc",times=dim(x)[1]))
soil[x[,"entvert"]==1] <- "entvert"        # 'Entisols','Vertisols'
soil[x[,"mol"]==1] <- "mol"                # 'Mollisols’
soil[x[,"spodhist"]==1] <- "spodhist"      # 'Spodosols','Histosols'
soil[x[,"ult"]==1] <- "ult"                # 'Ultisols_Udults_Others','Ultisols_Aquults','Ultisols_Others'
soil[x[,"ultkan"]==1] <- "ultkan"          # 'Ultisols_Udults_Kanhapludults'
soil[x[,"others"]==1] <- "others"          # 

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


# All
#X <- model.matrix(~ bs(stdage, 5) + bs(temp, 5) + bs(deficit, 5) + 
#                    bs(topo1, 5) + bs(topo2, 5) + bs(topo3, 5) + as.factor(soil), 
#                  data=rbind(dat.train, dat.test))
#Y <- c(dat.train$baBlackCherry, dat.test$baBlackCherry)==0

# Lasso
# Cross-validation on all data to pick lambda
# fit.ls <- cv.glmnet(x=X, y=Y, family="binomial")  

# Train and test
# fit.ls <- glmnet(x=X.train, y=Y.train, family="binomial", lambda=fit.ls$lambda.min)
# auc(Y.test, as.vector(predict(fit.ls, X.test, type="response")))

#####################################################################################
##################################  MODEL FIT  ######################################
#####################################################################################
# Functions for metric calculations
  
  source("~/Dropbox/Duke-Research/sims/zeroinf/0ComparisonFunctions.R")  

##########
# Function for running Cross Validated AUC metric for Elastic Net Model
do.loss.en.work <- function(j, data.title, fit.title){
  
  load(file = paste("kfold/", data.title, formatC(j, width = 2, format = "d", flag = "0"), ".Rdata", sep = ""))
  
  fac.soil.train <- model.matrix(~-1 + as.factor(soil), data=dat.train)
  fac.soil.test <- model.matrix(~-1 + as.factor(soil), data=dat.test)
  bs.stdage.train <- bs(dat.train$stdage, 5)
  bs.stdage.test <- predict(bs.stdage.train, dat.test$stdage)
  bs.temp.train <- bs(dat.train$temp, 5)
  bs.temp.test <- predict(bs.temp.train, dat.test$temp)
  bs.deficit.train <- bs(dat.train$deficit, 5)
  bs.deficit.test <- predict(bs.deficit.train, dat.test$deficit)
  bs.topo1.train <- bs(dat.train$topo1, 5)
  bs.topo1.test <- predict(bs.topo1.train, dat.test$topo1)
  bs.topo2.train <- bs(dat.train$topo2, 5)
  bs.topo2.test <- predict(bs.topo2.train, dat.test$topo2)
  bs.topo3.train <- bs(dat.train$topo3, 5)
  bs.topo3.test <- predict(bs.topo3.train, dat.test$topo3) 
  
  # Train
  X.train <- cbind(bs.stdage.train, bs.temp.train, bs.deficit.train,
                   bs.topo1.train, bs.topo2.train, bs.topo3.train, fac.soil.train)
  Y.train <- (dat.train$baBlackCherry==0)

  # Test
  X.test <- cbind(bs.stdage.test, bs.temp.test, bs.deficit.test,
                  bs.topo1.test, bs.topo2.test, bs.topo3.test, fac.soil.test)
  Y.test <- (dat.test$baBlackCherry==0)
  
    
  ###########################################################################################
  # BINOMIAL ELASTIC NET
  
  # Fixed for model with linear-linear splines on topo1 and temp
  a <- seq(0.1, 0.9, 0.05)
  cvm <- lambda.1se <- rep(NA, length(a)) 
  for (i in 1:length(a)) {
    cv <- cv.glmnet(X.train, Y.train, family = "binomial", nfold = 10, type.measure = "deviance", alpha = a[i])
    cvm[i] <- cv$cvm[cv$lambda == cv$lambda.1se]
    lambda.1se[i] <- cv$lambda.1se
  }
  lambda.use <- lambda.1se[which.min(cvm)]
  alpha.use <- a[which.min(cvm)]

  # Train and test
  fit.en <- glmnet(X.train, Y.train, family = "binomial", lambda = lambda.use, alpha = alpha.use)
  prob0.en.test <- as.vector(predict(fit.en, newx=X.test, type="response"))
  
  ###########################################################################################
  # COMPARISON CALCULATIONS
  
  prob0 <- data.frame(auc.lg.test = auc(Y.test, prob0.en.test),
                      cre.lg.test = cross.ent(Y.test, prob0.en.test))
  
  write.table(prob0,file = paste("kfold/results_cens", formatC(j, width = 2, format = "d", flag = "0"), "_prob0glmnet",fit.title,".txt", sep = ""))
} 

mclapply(1:10, do.loss.en.work, data.title="paper", fit.title="paper_linlin")


##########  READ IN AUC RESULTS
  
  # Read in results of all simulations
  nsim <- 10
  fit.title <- "paper_linlin"
  readin <- function(jj) cbind(sim=jj,
                read.table(file = paste(loc, "kfold/results_cens", formatC(jj, width = 2, format = "d", flag = "0"), "_prob0glmnet", fit.title,".txt", sep = "")))
  sim_prob0 <- do.call("rbind",lapply(1:nsim,readin))
  prob0 <- matrix(colMeans(sim_prob0[,-1],na.rm=T),ncol=2,byrow=T)
  dimnames(prob0) <- list(c("en"), c("auc", "cre"))
  prob0
  #auc       cre
  #en 0.8163105 0.3625057

