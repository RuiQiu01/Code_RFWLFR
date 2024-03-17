##Frechet regression for spd##############

library("frechet")
library("shapes")
library("stringr")
library("pbapply")
library("parallel")
library("MASS")
library("matrixcalc")
source("LocCovReg_debug.R")
source("shape_revise.R")
source("FRFPackage2.R")
source("main.R")

set.seed(1)
run <- 100  #number of main loops
ISE <- matrix(0,run,4)

##Model settings################################################
setting <- 2
m <- 3  #dimension of covariance matrices
p <- 2  #number of features
n <- 100  #number of training data
n_xout <- 100  #number of testing data
method <- "Log_cholesky"; metric <- "log_cholesky"; alpha <- 1 #type of metric
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
  #the generation of Y
  M_obs=array(0,c(m,m,n))
  if (setting ==1) {
    for(j in 1:n){
      D <- matrix(c(1,cos((as.numeric(t(beta)%*%X_obs[j,]))*4*pi),cos((as.numeric(t(beta)%*%X_obs[j,]))*4*pi),1),2,2)
      Z <- matrix(0,2,2)
      Z[1,1] <- rnorm(1,0,1)
      Z[2,2] <- rnorm(1,0,1)
      Z[1,2]=Z[2,1] <- rnorm(1,0,1/sqrt(2))
      logY <- 0.2*Z+D
      M_obs[,,j] <-matrix.exp(logY) 
    }
  }
  if (setting ==2) {
    for(j in 1:n){
      pho_1 <- 0.8*cos(as.numeric(t(beta_1)%*%X_obs[j,])*4*pi)
      pho_2 <- 0.4*cos(as.numeric(t(beta_2)%*%X_obs[j,])*4*pi)
      D <- matrix(c(1,pho_1,pho_2,pho_1,1,pho_1,pho_2,pho_1,1),3,3)
      Z <- matrix(0,3,3)
      Z[1,1] <- rnorm(1,0,1)
      Z[2,2] <- rnorm(1,0,1)
      Z[3,3] <- rnorm(1,0,1)
      Z[1,2]=Z[2,1] <- rnorm(1,0,1/sqrt(2))
      Z[1,3]=Z[3,1] <- rnorm(1,0,1/sqrt(2))
      Z[2,3]=Z[3,2] <- rnorm(1,0,1/sqrt(2))
      logY <- 0.2*Z+D
      M_obs[,,j] <-matrix.exp(logY) 
    }
  }
  
  #GFR
  res1 <- GloCovReg(x=X_obs,M=M_obs,xout=xout,optns=list(corrOut=FALSE,metric=metric,alpha=alpha))
  #LFR(may fail)
  res2 <- NA
  if (p<=2){
    try(
      {res2 <- LocCovReg(x=X_obs,M=M_obs,xout=xout,optns=list(corrOut=FALSE,metric=metric,alpha=alpha))}
    )
  }
  #random forest weighted local frechet regression#####################
  #preparation of data format
  X <- list()
  X$type <- "scalar"
  X$id <- 1:n
  X$X <- X_obs
  Y <- list()
  Y$type <- "covariance"
  Y$id <- 1:n
  Y$Y <- M_obs
  Xout <- list()
  Xout$type <- "scalar"
  Xout$id <- 1:n_xout
  Xout$X <- xout
  #parallelized Frechet trees
  rf_shape_para <- function(Curve=NULL, Scalar=NULL, Factor=NULL,Shape=NULL,Image=NULL,Y,mtry,ntree, ncores,ERT=FALSE, aligned.shape=FALSE,ntry=3,timeScale=0.1, ...){
    
    cl <- makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    clusterExport(cl, c("X", "Y","deep","mtry","method","alpha", "bb","ntree"))
    clusterEvalQ(cl, library("stringr"))
    clusterEvalQ(cl, library("shapes"))
    clusterEvalQ(cl, source("shape_revise.R"))
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
  M_true=array(0,c(m,m,n_xout))
  if (setting ==1) {
    for(i in 1:n_xout){
      M_true[,,i] <- matrix.exp(matrix(c(1,cos((as.numeric(t(beta)%*%xout[i,]))*4*pi),cos((as.numeric(t(beta)%*%xout[i,]))*4*pi),1),2,2))
    }
  }
  if (setting ==2) {
    for(i in 1:n_xout){
      pho_1 <- 0.8*cos(as.numeric(t(beta_1)%*%xout[i,])*4*pi)
      pho_2 <- 0.4*cos(as.numeric(t(beta_2)%*%xout[i,])*4*pi)
      D <- matrix(c(1,pho_1,pho_2,pho_1,1,pho_1,pho_2,pho_1,1),3,3)
      M_true[,,i] <- matrix.exp(D)
    }
  }
  dist1=dist2=dist3=dist4 <- rep(0,n_xout)
  #test errors of GFR, RFWLCFR, RFWLLFR
  for (i in 1:n_xout){
    dist1[i] <- distcov(res1$Mout[[i]], M_true[,,i], method=method, alpha=alpha)^2
    dist3[i] <- distcov(res3$res[,,i], M_true[,,i], method=method, alpha=alpha)^2
    dist4[i] <- distcov(res4$res[,,i], M_true[,,i], method=method, alpha=alpha)^2
  }
  ISE1 <- mean(dist1) 
  ISE3 <- mean(dist3)
  ISE4 <- mean(dist4) 
  #test error of LFR
  ISE2 <- NA
  if (length(res2)>1){
    dist2 <- rep(0,n_xout)
    for (i in 1:n_xout){
      dist2[i] <- distcov(res2$Mout[[i]], M_true[,,i], method=method, alpha=alpha)^2
    }
    ISE2 <- mean(dist2) #calculate MSE
  }
  #record the results of 100 runs
  ISE[bb,] <- c(ISE1,ISE2,ISE3,ISE4)
}

MISE <- colMeans(ISE,na.rm = TRUE)
seISE <- apply(ISE, 2, sd, na.rm = TRUE)





