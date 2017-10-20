##########
# Function to run joint quantile regression. Needs input "dat.train" found in .Rdata file
do.work <- function(j, data.title, fit.title, design, nburn=1000, parstart="prior"){
  load(file = paste("kfold/", data.title, formatC(j, width = 2, format = "d", flag = "0"), ".Rdata", sep = ""))
  fit.gp <- qrjoint(design, data=dat.train, cens=ifelse(dat.train[,as.character(design[[2]])]==0,2,0),
                    nsamp=nburn, thin=20, par=parstart)  # burn-in
  fit.gp <- update(fit.gp, nadd=500, append=FALSE) # retain
  
  save(fit.gp, file = paste("kfold/fit", formatC(j, width = 2, format = "d", flag = "0"), fit.title, ".Rdata", sep = ""))
}

# Function to update joint quantile regression.
do.update <- function(j, fit.title, addn=500){
  load( file = paste("kfold/fit", formatC(j, width = 2, format = "d", flag = "0"), fit.title, ".Rdata", sep = ""))
  fit.gp <- update(fit.gp, nadd=addn, append=FALSE) # retain only latest adds
  save(fit.gp, file = paste("kfold/fit", formatC(j, width = 2, format = "d", flag = "0"), fit.title, ".Rdata", sep = ""))
}  

##########
# Function for running check loss comparison
do.loss.work <- function(j, data.title, fit.title, design){
  
  load(file = paste("kfold/", data.title, formatC(j, width = 2, format = "d", flag = "0"), ".Rdata", sep = ""))
  
  tau <- round(seq(0.01, 0.99, by=0.01),2)
  tau_ <- round(seq(0,1,by=0.01), 2)
  X.test  <- model.matrix(design,data=dat.test)
  p <- ncol(X.test)
  resp <- as.character(design[[2]])
  
  environment(design) <- environment()
  term <- attributes(terms(design))$term.labels
  rhs <- paste(attributes(terms(design))$term.labels, collapse=" + ")
  which.factors <- union(grep("as.factor",term),grep("poly",term)) # NOTE: if poly specified... logistic will do poly over bs
  if(length(which.factors)>0 & (length(term) - length(which.factors))>0) {
    rhs.lg <- paste(c(paste("bs(",term[-which.factors],",5)", sep=""), term[which.factors]), collapse=" + ")
  } else {if(length(which.factors)>0) {
    rhs.lg <- paste(term, collapse=" + ")
  } else {rhs.lg <- paste(paste("bs(",term,",5)", sep=""), collapse=" + ")}
  }
  design.lg <- formula(paste(design[[2]],"==0 ~ ", rhs.lg, sep=""))
  
  ###########################################################################################
  # JOINT QUANTILE REGRESSION
  
  load(file = paste("kfold/fit", formatC(j, width = 2, format = "d", flag = "0"), fit.title,".Rdata", sep = ""))
  beta.hat <- coef(fit.gp, plot = FALSE, burn=0, nmc = 500)   
  beta.gp <- beta.hat$beta.est[2:100,,]
  
  # Get quantile estimates
  pred.gp.test <- X.test %*% t(beta.gp[,,2])  # Piggy-back off of work already done by coef.
  pred.gp.test.0  <- pmax(pred.gp.test,0)
  dimnames(pred.gp.test.0) <- NULL
  prob0.gp.test <- apply(pred.gp.test, 1, function(f) approx(f, tau, 0, rule=2)$y)
  
  # Get quantile estimates when conditioned on Y > 0 
  Q.gp <- X.test %*%t(beta.hat$beta.est[,,2])
  Q.gp.cens <- array(NA,dim(Q.gp))
  for (i in 1:dim(X.test)[1]){
    tau_new <- tau_*(1 - prob0.gp.test[i]) + prob0.gp.test[i]
    Q.gp.cens[i,] <- unlist(lapply(tau_new, function(f) approx(tau_, Q.gp[i,], f, rule=1)$y))
  }
  Q.gp.cens <- Q.gp.cens[,-c(1,101)]
  
  ###########################################################################################
  # LOGISTIC REGRESSION
  
  # Training model
  fit.lg <- glm(design.lg, data = dat.train, family="binomial")
  prob0.lg.test <- predict(fit.lg, newdata=dat.test, type="response")
  
  ###########################################################################################
  # TOBIT REGRESSION
  
  # Training model
  fit.tb <- tobit(design, data = dat.train)
  
  # TOBIT Predictions: means dictate the entire functional curve
  pred.tb.test  <- sapply(tau, Vectorize(function(f) qnorm(f, predict(fit.tb, newdata=dat.test), fit.tb$scale)))
  pred.tb.test.0  <- pmax(pred.tb.test, 0)
  prob0.tb.test  <- pnorm(0, predict(fit.tb, newdata=dat.test), fit.tb$scale)
  
  # Conditional on Y>0
  Q.tb.cens <- array(NA,dim(pred.tb.test))
  for (i in 1:dim(X.test)[1]){
    tau_new <- tau*(1 - prob0.tb.test[i]) + prob0.tb.test[i]
    Q.tb.cens[i,] <- qnorm(tau_new, predict(fit.tb, newdata=dat.test[i,]), fit.tb$scale)
  }
  
  ###########################################################################################
  # ZERO-INFLATED (CENSORED) QUANTILE REGRESSION
  
  response.train <- survival::Surv(dat.train[,resp], dat.train[,resp]!=0, type = "left") 
  fit.pt <- crq(as.formula(paste("response.train ~", rhs, sep="")), data = dat.train, method="Portnoy")
  Sfit <- summary(fit.pt, tau)  # Does some sort of bootstrap or jackknife (not original)
  tau.pt <- round(sapply(Sfit, function(lo) return(lo$tau)),2)
  
  # Get quantile estimates
  betas.pt <- sapply(Sfit, function(lo) return(lo$coefficients[,1]))  # on original grid for obtaining quantiles
  pred.pt.test  <- X.test  %*% betas.pt
  pred.pt.test.0  <- pmax(pred.pt.test,0)
  
  stepfun.test  <- apply(pred.pt.test, 1, function(f){stepfun(tau.pt, c(f[1],f))})
  ra.test <- rearrange(stepfun.test)
  pred.pt.test.ra <- do.call("rbind", lapply(ra.test, function(f){environment(f)$y}))
  pred.pt.test.ra  <- pred.pt.test.ra[,-c(1,ncol(pred.pt.test.ra))]
  
  pred.pt.test.ra.0  <- pmax(pred.pt.test.ra,0)
  prob0.pt.test.ra <- apply(pred.pt.test.ra, 1, function(f) approx(f, tau.pt, 0, rule=1:2)$y)
  
  # Conditioning on Y >0
  prob0.pt.test.ra_ <- prob0.pt.test.ra
  prob0.pt.test.ra_[is.na(prob0.pt.test.ra_)] <- runif(sum(is.na(prob0.pt.test.ra_)), 0, min(tau.pt))
  Q.pt.cens <- array(NA,dim(pred.pt.test.ra))
  for (i in 1:dim(X.test)[1]){
    tau_new <- tau.pt*(1 - prob0.pt.test.ra_[i]) + prob0.pt.test.ra_[i]
    Q.pt.cens[i,] <- unlist(lapply(tau_new, function(f) approx(tau.pt, pred.pt.test.ra[i,], f, rule=1:2)$y)) #NA and closest for upper tail
  }
  
  # put portnoy estimates back on full tau grid
  padit <- length(tau) - length(which(tau %in% tau.pt))
  pred.pt.test <- cbind(matrix(NA,nrow(pred.pt.test),padit), pred.pt.test)
  pred.pt.test.0 <- cbind(matrix(NA,nrow(pred.pt.test.0),padit), pred.pt.test.0)
  pred.pt.test.ra <- cbind(matrix(NA,nrow(pred.pt.test.ra),padit), pred.pt.test.ra)
  pred.pt.test.ra.0 <- cbind(matrix(NA,nrow(pred.pt.test.ra.0),padit), pred.pt.test.ra.0)
  Q.pt.cens <- cbind(matrix(NA,nrow(Q.pt.cens),padit), Q.pt.cens)
  
  
  ###########################################################################################
  ###########################################################################################
  # CHECK LOSS CALCULATIONS
  
  cl <- data.frame(tau,
                   Tobit=colMeans(check.loss(dat.test[,resp] - pred.tb.test.0, tau)),
                   Portnoy=colMeans(check.loss(dat.test[,resp] - pred.pt.test.0, tau)),
                   Chernoz=colMeans(check.loss(dat.test[,resp] - pred.pt.test.ra.0, tau)),
                   JointQR=colMeans(check.loss(dat.test[,resp] - pred.gp.test.0, tau)))
  
  cl_cens <- data.frame(tau,
                        Tobit=colMeans(check.loss((dat.test[,resp] - Q.tb.cens)[dat.test[,resp]>0,], tau)),
                        Portnoy=colMeans(check.loss((dat.test[,resp] - Q.pt.cens)[dat.test[,resp]>0,], tau)),
                        JointQR=colMeans(check.loss((dat.test[,resp] - Q.gp.cens)[dat.test[,resp]>0,], tau)))
  
  prob0 <- data.frame(auc.lg.test = auc(dat.test[,resp]==0, prob0.lg.test),
                      cre.lg.test = cross.ent(dat.test[,resp]==0, prob0.lg.test),
                      auc.tb.test = auc(dat.test[,resp]==0, prob0.tb.test),
                      cre.tb.test = cross.ent(dat.test[,resp]==0, prob0.tb.test),
                      auc.pt.test = auc.mc(dat.test[,resp]==0, prob0.pt.test.ra, tau=tau),
                      cre.pt.test = cross.ent.mc(dat.test[,resp]==0, prob0.pt.test.ra, tau=tau),
                      auc.gp.test = auc(dat.test[,resp]==0, prob0.gp.test),
                      cre.gp.test = cross.ent(dat.test[,resp]==0, prob0.gp.test))
  
  write.table(cl,file = paste("kfold/results_cens", formatC(j, width = 2, format = "d", flag = "0"), "_cl",fit.title,".txt", sep = ""))
  write.table(cl_cens, file = paste("kfold/results_cens", formatC(j, width = 2, format = "d", flag = "0"), "_cl_cens",fit.title,".txt", sep = ""))
  write.table(prob0,file = paste("kfold/results_cens", formatC(j, width = 2, format = "d", flag = "0"), "_prob0",fit.title,".txt", sep = ""))
} 


##########
# Function to plot results. Will only work if it finds appropriately named files.
result.plots <- function(nsim, nsub=NULL, fit.title){
  
  if(!is.null(nsim)){
    use <- 1:nsim
  }
  if(!is.null(nsub)){
    use <- nsub
  }
  # Read in results of all simulations
  readin <- function(jj,summ) cbind(sim=jj,
                                    read.table(file = paste(loc, "kfold/results_cens", formatC(jj, width = 2, format = "d", flag = "0"), summ, sep = "")))
  
  sim_cl <- do.call("rbind",lapply(use,readin, summ=paste("_cl", fit.title,".txt",sep="")))
  sim_cl_cens <- do.call("rbind",lapply(use,readin, summ=paste("_cl_cens", fit.title,".txt",sep="")))
  sim_prob0 <- do.call("rbind",lapply(use,readin,summ=paste("_prob0", fit.title,".txt",sep="")))
  
  #####################
  # Perform Calculations
  
  # Calculations for relative efficiency check loss (censored version)
  cl <- sim_cl %>% 
    group_by(tau) %>% 
    summarize(Tobit=mean(Tobit),Portnoy=mean(Portnoy),
              JointQR=mean(JointQR)) %>% as.data.frame
  cl_rel <- cbind(tau=cl$tau,
                  JointQR = cl$Tobit/cl$JointQR,
                  Portnoy = cl$Tobit/cl$Portnoy)
  cl_rel <- as.data.frame(cl_rel) %>% gather(method,CheckLoss, -tau)
  
  # Calculations for relative efficiency check loss (Y>0 conditioning)
  cl_cens <- sim_cl_cens %>% 
    group_by(tau) %>% 
    summarize(Tobit=mean(Tobit),
              Portnoy=mean(Portnoy),
              JointQR=mean(JointQR)) %>% as.data.frame
  cl_cens_rel <- cbind(tau=cl_cens$tau,
                       JointQR = cl_cens$Tobit/cl_cens$JointQR,
                       Portnoy = cl_cens$Tobit/cl_cens$Portnoy)
  cl_cens_rel <- as.data.frame(cl_cens_rel) %>% gather(method,CheckLoss, -tau)
  
  # Calculations for probability of 0
  prob0 <- matrix(colMeans(sim_prob0[,-1],na.rm=T),ncol=2,byrow=T)
  dimnames(prob0) <- list(c("lg","tb","pt","gp"), c("auc", "cre"))
  
  #####################
  # Plot results
  
  # Relative Censored Checkloss Compared to Tobit
  p.cens <- ggplot(data=cl_rel, aes(x=tau, y=CheckLoss, col=method, group=method)) + geom_line() +
    geom_abline(intercept=1, slope=0, col="gray") + 
    xlab(expression(paste("Response proportion, ",tau))) + ylab("Efficency Relative to Tobit") + 
    scale_color_manual(values=c("blue", "red")) + ggtitle("Censored Check Loss")
  
  # Relative Checkloss Compared to Tobit (Truncated Quantiles)
  p.trunc <- ggplot(data=cl_cens_rel, aes(x=tau, y=CheckLoss, col=method, group=method)) + geom_line() +
    geom_abline(intercept=1, slope=0, col="gray") + 
    xlab(expression(paste("Response proportion, ",tau))) + ylab("Efficency Relative to Tobit") +
    scale_color_manual(values=c("blue", "red")) + ggtitle("Truncated Quantile Check Loss")
  
  pdf(file=paste("~/Dropbox/Duke-Research/casestudies/environment/CrossValCheckLoss",fit.title,".pdf",sep=""), width=8, height=3)
  grid.arrange(p.cens, p.trunc, ncol=2) 
  dev.off()
  
  prob0
}

