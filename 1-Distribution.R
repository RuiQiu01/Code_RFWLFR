##Frechet regression for distributions##############

library("frechet")
library("stringr")
library("pbapply")
library("parallel")
library("MASS")
source("FRFPackage2.R")
source("main.R")

set.seed(1)
run <- 100  #number of main loops
ISE <- matrix(0,run,4)

##Model settings################################################
setting <- 2
p <- 2   #dimension of features
n <- 100   #number of training data
n_xout <- 100  #number of testing data
if (setting ==1) {
  if (p==2) {beta <- c(0.75,0.25)}
  if (p==5) {beta <- c(0.1,0.2,0.3,0.4,0)}
  if (p==10) {beta <- c(0.1,0.2,0.3,0.4,rep(0,6))}
  if (p==20) {beta <- c(0.1,0.2,0.3,0.4,rep(0,12),0.1,0.2,0.3,0.4)/2}
}
if (setting ==2) {
  if (p==2) {beta_1 <- c(0.75,0.25);beta_2 <- c(0.25,0.75)}
  if (p==5) {beta_1 <- c(0.1,0.2,0.3,0.4,0);beta_2 <- c(0,0.1,0.2,0.3,0.4)}
  if (p==10) {beta_1 <- c(0.1,0.2,0.3,0.4,rep(0,6));beta_2 <- c(rep(0,6),0.1,0.2,0.3,0.4)}
  if (p==20) {beta_1 <- c(0.1,0.2,0.3,0.4,rep(0,16));beta_2 <- c(rep(0,16),0.1,0.2,0.3,0.4)}
}

##main loops####################################################
for (bb in 1:run) {
  data <- runif(n*p, 0, 1) 
  X_obs <- matrix(data, n, p) #trainning data
  data_out <- runif(n_xout*p, 0, 1)
  xout <- matrix(data_out , n_xout, p) #testing data
  nqSup <- 21 
  qSup <- seq(0,1,length.out = nqSup)
  #the generation of Y
  mu <- rep(0,n)
  sigma <- rep(0,n)
  Y_obs<- matrix(0,n,nqSup) 
  if (setting ==1) {
    for (i in 1:n) {
      a <- as.numeric(t(beta)%*%X_obs[i,])
      mu[i] <- rnorm(1,5*a-2.5,0.2)
      sigma[i] <- 1
      Y_obs[i,] <- qnorm(qSup,mu[i],sigma[i])
    }
    Y_obs[,1] <- Y_obs[,2]
    Y_obs[,nqSup] <- Y_obs[,nqSup-1]  
  }
  if (setting ==2) {
    for (i in 1:n) {
      a <- as.numeric(t(beta_1)%*%X_obs[i,]);b <- as.numeric(t(beta_2)%*%X_obs[i,])
      mu[i] <- rnorm(1,sin(4*pi*a)*(2*b-1),0.2)
      sigma[i] <- 2*abs(X_obs[i,1]-X_obs[i,2])
      Y_obs[i,] <- qnorm(qSup,mu[i],sigma[i])
    }
    Y_obs[,1] <- Y_obs[,2]
    Y_obs[,nqSup] <- Y_obs[,nqSup-1] 
  }
  
  #GFR
  res1 <- GloDenReg(xin=X_obs, qin=Y_obs, xout=xout, optns = list(qSup = qSup))
  #LFR(may fail)
  res2 <- NA
  if (p<=2){
    try(
      {res2 <- LocDenReg(xin=X_obs, qin=Y_obs, xout=xout, optns = list(qSup = qSup))}
    )
  }
  #random forest weighted local frechet regression#####################
  #preparation of data format
  X <- list()
  X$type <- "scalar"
  X$id <- 1:n
  X$X <- X_obs
  Y <- list()
  Y$type <- "distribution"
  Y$id <- 1:n
  Y$Y <- Y_obs
  Xout <- list()
  Xout$type <- "scalar"
  Xout$id <- 1:n_xout
  Xout$X <- xout
  #parallelized Frechet trees
  rf_shape_para <- function(Curve=NULL, Scalar=NULL, Factor=NULL,Shape=NULL,Image=NULL,Y,mtry,ntree, ncores,ERT=FALSE, aligned.shape=FALSE,ntry=3,timeScale=0.1, ...){

    cl <- makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    clusterExport(cl, c("X", "Y","qSup","deep","mtry","nqSup","bb","ntree"))
    clusterEvalQ(cl, library("stringr"))
    clusterEvalQ(cl, source("FRFPackage2.R"))
    
    trees <- pbsapply(1:ntree, FUN=function(i){
      set.seed(i+bb*ntree)
      Rtmax(Curve=Curve,Scalar = Scalar,Factor = Factor,Shape=Shape,Image=Image,Y,mtry,ERT=ERT, aligned.shape=aligned.shape,ntry=ntry,timeScale=timeScale, ...)
    },cl=cl)
    
    stopCluster(cl)
    
    return(trees)
  }
  #RFWLCFE
  res3 <- rfwlcfr(Scalar=X, Y=Y, Xout=Xout, mtry=ceiling(4/5*p), deep=2, ntree=100, ncores=50)
  #RFWLLFR
  res4 <- rfwllfr(Scalar=X, Y=Y, Xout=Xout, mtry=ceiling(4/5*p), deep=5, ntree=100, ncores=50)
  
  #measure testing error#################################################
  #generate true responses for testing data
  mu_true <- rep(0,n_xout) 
  sigma_true <- rep(0,n_xout)
  Y_true<- matrix(0,n_xout,nqSup)
  if (setting ==1) {
    for (i in 1:n_xout) {
      a <- as.numeric(t(beta)%*%xout[i,])
      mu_true[i] <- 10*xout[i,1]^2*(2*a-1)
      sigma_true[i] <- 1
      Y_true[i,] <- qnorm(qSup, mu_true[i], sigma_true[i])
    }
    Y_true[,1] <- Y_true[,2]
    Y_true[,nqSup] <- Y_true[,nqSup-1]
  }
  if (setting ==2) {
    for (i in 1:n_xout) {
      a <- as.numeric(t(beta_1)%*%xout[i,]);b <- as.numeric(t(beta_2)%*%xout[i,])
      mu_true[i] <- sin(4*pi*a)*(2*b-1)
      sigma_true[i] <- 2*abs(xout[i,1]-xout[i,2])
      Y_true[i,] <- qnorm(qSup, mu_true[i], sigma_true[i])
    }
    Y_true[,1] <- Y_true[,2]
    Y_true[,nqSup] <- Y_true[,nqSup-1]
  }
  dist1=dist2=dist3=dist4 <- rep(0,n_xout)
  #test errors of GFR, RFWLCFR, RFWLLFR
  for (i in 1:n_xout){
    dist1[i] <- pracma::trapz(qSup, (res1$qout[i,] - Y_true[i,])^2)
    dist3[i] <- pracma::trapz(qSup, (res3$res[i,] - Y_true[i,])^2)
    dist4[i] <- pracma::trapz(qSup, (res4$res[i,] - Y_true[i,])^2)
  }
  ISE1 <- mean(dist1)
  ISE3 <- mean(dist3) 
  ISE4 <- mean(dist4)
  #test error of LFR
  ISE2 <- NA
  if (length(res2)>1){
    for (i in 1:n_xout){
      dist2[i] <- pracma::trapz(qSup, (res2$qout[i,] - Y_true[i,])^2)
    }
    ISE2 <- mean(dist2) 
  }
  #record the results of 100 runs
  ISE[bb,] <- c(ISE1,ISE2,ISE3,ISE4)
}

MISE <- colMeans(ISE,na.rm = TRUE)
seISE <- apply(ISE, 2, sd, na.rm = TRUE)




