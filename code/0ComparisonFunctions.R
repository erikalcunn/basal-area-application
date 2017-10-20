###########################################
# Functions common to multiple simulations or k-fold validation schemes


##########
# Function to calculate check.loss function
# takes n x L "res" predicted quantile residuals as input, and 
#              tau (of length L and coorsponding to dim of res) as input
check.loss <- function(res,tau){res*(matrix(tau, dim(res)[1],length(tau),byrow=T) - 1*(res<0))}







#########
# Function to estimate cross entropy
#  -sum(y_i *log(phat_i) + (1 - y_i)*log(1-phat_i))
# Take as input actual: binary 01 data
#                 pred: predicted probability of being zero
cross.ent <- function(actual,pred, tau=NULL){-1*mean(actual*log(pred) + (1-actual)*log(1-pred))}  

#########
# Function to estimate cross entropy in the case of deprecation of the
# prediction values to  P(X==0) < min(tau). Indicated by NAs in the predictions
# Uses 100 monte carlo draws to do so
#     tau:  (for deprecated prob 0 predictions, a single value or vector
#           over which probabilities were able to be estimated
#           when tau given -- predicted probabilities are sampled uniformly from 0 to min(tau)
cross.ent.mc <- function(actual,pred, tau=NULL){
  cross.ent.draw <- function(actual.=actual, pred.=pred, tau.=tau){
    pred.temp <- pred.
    pred.temp[is.na(pred)] <- runif(sum(is.na(pred.)), 0, min(tau.))
    -1*mean(actual*log(pred.temp) + (1-actual)*log(1-pred.temp))
  }
  mean(sapply(1:100, cross.ent.draw))
} 

##########
# A function that uses 100 monte carlo draws to estimate AUC when the predictions for
# zero-inflaction are deprecated such that P(X==0) < min(tau) but otherwise unknown
# NAs in the pred vector indicate the deprecation.
auc.mc <- function(actual, pred, tau){
  auc.mc.draw <- function(actual.=actual, pred.=pred, tau.=tau){
    pred.temp <- pred.
    pred.temp[is.na(pred)] <- runif(sum(is.na(pred.)), 0, min(tau.))
    auc(actual, pred.temp)
  }
  mean(sapply(1:100, auc.mc.draw))
}

##########
# Function for evaluating sensitivity, specificity, accuracy
classperf <- function(actual, prob.0, tau, plot=TRUE){
  sens <- spec <- accr <- rep(NA,length(tau))
  for (i in tau){
    tab <- table(actual, prob.0> i)
    if(all(dim(tab)==c(2,2))){
      
      sens[which(tau==i)] <- prop.table(tab,1)[2,2]
      spec[which(tau==i)] <- prop.table(tab,1)[1,1]
      accr[which(tau==i)] <- (tab[1,1] + tab[2,2])/sum(tab)
    }
  }
  if(plot){
    temp.par <- par()$mfrow
    par(mfrow=c(2,2))
    plot(tau, sens, type="s", xlim=c(0,1), ylim=c(0,1))
    plot(tau, spec, type="s", xlim=c(0,1), ylim=c(0,1))
    plot(tau, accr, type="s", xlim=c(0,1), ylim=c(0,1))
    plot(1-spec, sens, type="s", xlim=c(0,1), ylim=c(0,1))
    abline(0,1, col="gray")
    par(mfrow = temp.par)
  } else{
    return(list(sens=sens, spec=spec, accr=accr))}
}