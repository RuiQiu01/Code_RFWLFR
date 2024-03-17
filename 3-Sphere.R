##Frechet regression for spherical data##############

library("stringr")
library("pbapply")
library("parallel")
library("RiemBase")
library("spherepc")
library("MASS")
source("IntrinsicMean_revise.R")
source("FRFPackage2.R")
source("main.R")

set.seed(1)
run <- 100
ISE <- matrix(0,run,2)

##Model settings################################################
setting <- 1
p <- 2  #number of features
n <- 100  #number of training data
n_xout <- 100  #number of testing data
if (p==2) {beta_1 <- c(1,0);beta_2 <- c(0,1)}
if (p==5) {beta_1 <- c(0.1,0.2,0.3,0.4,0);beta_2 <- c(0,0.1,0.2,0.3,0.4)}
if (p==10) {beta_1 <- c(0.1,0.2,0.3,0.4,rep(0,6));beta_2 <- c(rep(0,6),0.1,0.2,0.3,0.4)}
if (p==20) {beta_1 <- c(0.1,0.2,0.3,0.4,rep(0,16));beta_2 <- c(rep(0,16),0.1,0.2,0.3,0.4)}

##main loops####################################################
bb <- 1
while (bb <= run) {
  fail <- 0
  data <- runif(n*p, 0, 1) 
  X_obs <- matrix(data, n, p) #trainning data
  data_out <- runif(n_xout*p, 0, 1)
  xout <- matrix(data_out , n_xout, p) #testing data
  #the generation of Y
  if (setting ==1) {
    S_t <- matrix(0,n,3)
    S_obs <- matrix(0,n,3)
    for (i in 1:n){
      tmpx <- sqrt((1-(t(beta_1)%*%X_obs[i,])^2))*cos(pi*(t(beta_2)%*%X_obs[i,]))
      tmpy <- sqrt((1-(t(beta_1)%*%X_obs[i,])^2))*sin(pi*(t(beta_2)%*%X_obs[i,]))
      tmpz <- t(beta_1)%*%X_obs[i,]
      S_t[i,] <- c(tmpx,tmpy,tmpz)
      A <- svd(diag(3)-S_t[i,]%*%t(S_t[i,]),nu=2)$u
      z <-  0.2*rnorm(2)
      epsilon <- A%*%z
      norm_epsilon <- sqrt(sum(epsilon^2))
      S_obs[i,] <- cos(norm_epsilon)*S_t[i,]+sin(norm_epsilon)*epsilon/norm_epsilon
    }
  }
  if (setting ==2) {
    S_obs <- matrix(0,n,3)
    for (i in 1:n){
      eps <- rnorm(2,0,0.2)
      tmpx <- sin(t(beta_1)%*%X_obs[i,]+eps[1])*sin(t(beta_2)%*%X_obs[i,]+eps[2])
      tmpy <- sin(t(beta_1)%*%X_obs[i,]+eps[1])*cos(t(beta_2)%*%X_obs[i,]+eps[2])
      tmpz <- cos(t(beta_1)%*%X_obs[i,]+eps[1])
      S_obs[i,] <- c(tmpx,tmpy,tmpz)
    }
  }
  
  #random forest weighted local frechet regression#####################
  #preparation of data format
  X <- list()
  X$type <- "scalar"
  X$id <- 1:n
  X$X <- X_obs
  Y <- list()
  Y$type <- "sphere"
  Y$id <- 1:n
  Y$Y <- S_obs
  Xout <- list()
  Xout$type <- "scalar"
  Xout$id <- 1:n_xout
  Xout$X <- xout
  #parallelized Frechet trees
  rf_shape_para <- function(Curve=NULL, Scalar=NULL, Factor=NULL,Shape=NULL,Image=NULL,Y,mtry,ntree, ncores,ERT=FALSE, aligned.shape=FALSE,ntry=3,timeScale=0.1, ...){
    
    cl <- makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    clusterExport(cl, c("X", "Y","deep","mtry","bb","ntree"))
    clusterEvalQ(cl, library("stringr"))
    clusterEvalQ(cl, library("RiemBase"))
    clusterEvalQ(cl, source("FRFPackage2.R"))
    
    trees <- pbsapply(1:ntree, FUN=function(i){
      set.seed(i+bb*ntree)
      Rtmax(Curve=Curve,Scalar = Scalar,Factor = Factor,Shape=Shape,Image=Image,Y,mtry,ERT=ERT, aligned.shape=aligned.shape,ntry=ntry,timeScale=timeScale, ...)
    },cl=cl)
    
    stopCluster(cl)
    
    return(trees)
  }
  #RFWLCFE
  res3 <- rfwlcfr(Scalar=X, Y=Y, Xout=Xout, mtry=ceiling(1/5*p), deep=3, ntree=100, ncores=50)
  if (is.numeric(res3)){next}
  #RFWLLFR
  res4 <- rfwllfr(Scalar=X, Y=Y, Xout=Xout, mtry=ceiling(1/5*p), deep=15, ntree=100, ncores=50)
  if (is.numeric(res3)){next}
  
  #measure testing error#################################################
  #generate true responses for testing data
  S_true=matrix(0,n_xout,3)
  if (setting ==1) {
    for(i in 1:n_xout){
      tmpx <- sqrt((1-(t(beta_1)%*%xout[i,])^2))*cos(pi*(t(beta_2)%*%xout[i,]))
      tmpy <- sqrt((1-(t(beta_1)%*%xout[i,])^2))*sin(pi*(t(beta_2)%*%xout[i,]))
      tmpz <- t(beta_1)%*%xout[i,]
      S_true[i,] <- c(tmpx,tmpy,tmpz)
    }
  }
  if (setting ==2) {
    for(i in 1:n_xout){
      tmpx <- sin(t(beta_1)%*%xout[i,])*sin(t(beta_2)%*%xout[i,])
      tmpy <- sin(t(beta_1)%*%xout[i,])*cos(t(beta_2)%*%xout[i,])
      tmpz <- cos(t(beta_1)%*%xout[i,])
      S_true[i,] <- c(tmpx,tmpy,tmpz)
    }
  }
  dist3=dist4 <- rep(0,n_xout)
  #test errors of RFWLCFR, RFWLLFR
  for (i in 1:n_xout){
    cross <- crossprod(res3$res[i,],S_true[i,])[1,1]
    dist3[i] <- acos(cross)^2
    cross <- crossprod(res4$res[i,],S_true[i,])[1,1]
    dist4[i] <- acos(cross)^2
  }
  ISE3 <- mean(dist3) 
  ISE4 <- mean(dist4) 
  #record the results of 100 runs
  ISE[bb,] <- c(ISE3,ISE4)
  bb <- bb+1
}

MISE <- colMeans(ISE,na.rm = TRUE)
seISE <- apply(ISE, 2, sd, na.rm = TRUE)


