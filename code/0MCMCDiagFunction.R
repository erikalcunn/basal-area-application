# Function for convergence checks ("psrf" potential scale reduction factor), modified from Surya,
# on untransformed parameters and on betas (finds max and multivariate for each).
do.diag <- function(j, fit.title, fit.title2, nsamp){
  
  # First run, plot diganostics, look at WAIC
  load(file = paste(loc,"kfold/fit", formatC(j, width = 2, format = "d", flag = "0"), fit.title, ".Rdata", sep = ""))
  o.list <- list(fit.gp)
  pdf(paste(loc,"kfold/diag", formatC(j, width = 2, format = "d", flag = "0"), fit.title, ".pdf", sep = ""), width=8, height=18)
  WAIC <- summary(fit.gp, more.details=TRUE)$waic
  dev.off()
  
  rm(fit.gp)
  
  load(file = paste("kfold/fit", formatC(j, width = 2, format = "d", flag = "0"), fit.title2, ".Rdata", sep = ""))
  o.list[[2]] <- fit.gp
  rm(fit.gp)
  
  # pass o.list, containing at least 2 qrjoint fitted objects
  pars.list <- mcmc.list(lapply(o.list, function(o) return(mcmc(t(matrix(o$parsamp, ncol = nsamp)))))) 
  
  # gelman diagnostics on untransformed parameters, two diffuse chains 
  pars.diag <- gelman.diag(pars.list)
  pars.psrf <- c(apply(pars.diag$psrf, 2, max), mpsrf = pars.diag$mpsrf)
  
  # gelman diagnostics on beta parameters, two chains.
  beta.samp.list <- mclapply(o.list, function(f) {coef(f, plot = FALSE, burn = 0, nmc = 500)$beta.samp})
  beta.list <- mcmc.list(lapply(beta.samp.list, function(bb) return(mcmc(t(rbind(bb[11,,],bb[51,,],bb[91,,]))))))
  beta.diag <- gelman.diag(beta.list)
  beta.psrf <- c(apply(beta.diag$psrf, 2, max), mpsrf = beta.diag$mpsrf) # First set 10%, next set 50%, Last set 90%
  return(list(beta.diag, pars.psrf, beta.psrf, WAIC))
}