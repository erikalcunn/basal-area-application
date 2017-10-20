##########  
# Function to split data into train/test sets and save data for later access.
# can subset data but will stratify within soil categories when doing the k-fold split

data.create <- function(nsim, subset, data.title){
  
  dat2 <- dat[subset, ]
  dat2$soil <- factor(dat2$soil)
  soil.n <- table(dat2$soil)
  p <- length(soil.n)
  
  # Full-length random sequences and repeated 1:nsim vectors to go along with each
  for (i in 1:p){
    assign(paste("ind",i,sep=""), sample(1:soil.n[i], soil.n[i], replace=FALSE))
    assign(paste("kfold",i,sep=""), rep(1:nsim, length=soil.n[i]))
  }
  
  for (j in 1:nsim){
    for (i in 1:p){
      assign(paste("train",i,sep=""), get(paste("ind",i,sep=""))[get(paste("kfold",i,sep=""))!=j])
      assign(paste("dat.train",i,sep=""), dat2[dat2$soil==names(soil.n[i]),][get(paste("train",i,sep="")),]) 
      assign(paste("dat.test",i,sep=""), dat2[dat2$soil==names(soil.n[i]),][-get(paste("train",i,sep="")),]) 
    }
    
    dat.train <- do.call("rbind",lapply(1:p, function(f) get(paste("dat.train",f,sep=""))))
    dat.test <- do.call("rbind",lapply(1:p, function(f) get(paste("dat.test",f,sep=""))))
    save(dat.train, dat.test, file = paste("kfold/", data.title, formatC(j, width = 2, format = "d", flag = "0"), ".Rdata", sep = ""))
  } 
}  