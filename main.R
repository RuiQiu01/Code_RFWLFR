#' @title Random forest weighted local constant Frechet regreession.
#' @param Scalar A list about input measurements of the predictors from training data.
#' @param Y A list about output measurements of the responses from training data, which can be distributions, covariances, spherical data.
#' @param Xout An n_out by p matrix with input measurements of the predictors from testing data.
#' @param mtry The number of feature subsampling during the tree building, a number between 1 and p.
#' @param deep The depth of tree growth. When using the package "FRFPackage", "deep" is the maximum number of split layers of the tree；When using the package "FRFPackage2", "deep" is the maximum number of samples within a leaf. Here the package "FRFPackage2" is used. 
#' @param ntree The number of Frechet trees in the forest. Default is 100.
#' @param ncores The number of cores used for parallel operations when building trees.
rfwlcfr <- function(Scalar=X, Y=Y, Xout=Xout, mtry=mtry, deep=2, ntree=100, ncores=ncores){
  mtry <<- mtry; deep <<- deep; ntree <<- ntree
  FRF <- FrechForest(Scalar=Scalar, Y=Y, mtry=mtry, ntree=ntree, ncores=ncores, ERT=FALSE, OOB=FALSE, imp=FALSE)
  res3 <- predict.FrechForest(object=FRF, Scalar=Xout)
  n_xout <- max(Xout$id); n <- max(X$id)
  RF_weight <- matrix(0,n_xout,n)
  for (i in 1:n_xout) {
    Tree_weight <- matrix(0,ntree,n)
    for (j in 1:ntree) {
      neighbor <- FRF$rf[,j]$idY[which(FRF$rf[,j]$feuilles==res3$pred_feuille[j,i])]
      for (k in 1:n) {
        if (k %in% neighbor) {
          Tree_weight[j,k] <- 1/length(neighbor)
        }
      }
    }
    RF_weight[i,] <-colMeans(Tree_weight) 
  }
  if (Y$type == "distribution"){
    res4 <- RF_weight %*% Y_obs}
  if (Y$type == "covariance"){
    res4 <- array(0,c(m,m,n_xout))
    for (i in 1:n_xout) {
      res4[,,i] <- estcov(M_obs, method=method, weights=RF_weight[i,], alpha=alpha)$mean
    }}
  fail <- 0
  if (Y$type == "sphere"){
    #converting 3d coordinates to latitude and longitude coordinates to calculate weighted mean
    S2_obs <- matrix(0,n,2)
    for (i in 1:n){
      S2_obs[i,] = Trans.sph(S_obs[i,]) 
    }
    Res4 <- matrix(0,n_xout,2)
    for (i in 1:n_xout) {
      Res4[i,] <- IntrinsicMean(S2_obs, weights = RF_weight[i,])
      if (Res4[i,1]==5) {fail <- 1;break}
    }
    if (fail==0) {
      #convert back to 3d coordinates
      res4 <- matrix(0,n_xout,3)
      for (i in 1:n_xout) {
        res4[i,] <- Trans.Euclid(Res4[i,])
      }
    }}
  if (fail==0) {return(list(res=res4, weight=RF_weight))}
  if (fail==1) {return(1)}
}

#' @title Random forest weighted local linear Frechet regreession.
#' @param Scalar An n by p matrix with input measurements of the predictors from training data.
#' @param Y Measurements of the corresponding responses from training data, which can be distributions, covariances, spherical data.
#' @param Xout An n_out by p matrix with input measurements of the predictors from testing data.
#' @param mtry The number of feature subsampling during the tree building, a number between 1 and p.
#' @param deep The depth of tree growth. When using the package "FRFPackage", "deep" is the maximum number of split layers of the tree；When using the package "FRFPackage2", "deep" is the maximum number of samples within a leaf. Here the package "FRFPackage2" is used. 
#' @param ntree The number of Frechet trees in the forest. Default is 100.
#' @param ncores The number of cores used for parallel operations when building trees.
rfwllfr <- function(Scalar=X, Y=Y, Xout=Xout, mtry=mtry, deep=5, ntree=100, ncores=ncores){
  mtry <<- mtry; deep <<- deep; ntree <<- ntree
  FRF <- FrechForest(Scalar=Scalar, Y=Y, mtry=mtry, ntree=ntree, ncores=ncores, ERT=FALSE, OOB=FALSE, imp=FALSE)
  res3 <- predict.FrechForest(object=FRF, Scalar=Xout)
  n_xout <- max(Xout$id); n <- max(X$id)
  RF_weight <- matrix(0,n_xout,n)
  for (i in 1:n_xout) {
    Tree_weight <- matrix(0,ntree,n)
    for (j in 1:ntree) {
      neighbor <- FRF$rf[,j]$idY[which(FRF$rf[,j]$feuilles==res3$pred_feuille[j,i])]
      for (k in 1:n) {
        if (k %in% neighbor) {
          Tree_weight[j,k] <- 1/length(neighbor)
        }
      }
    }
    RF_weight[i,] <-colMeans(Tree_weight) 
  }
  loc_weight <- matrix(0,n_xout,n)
  for (i in 1:n_xout) {
    D <- X_obs-matrix(rep(xout[i,],n),n,p,byrow = T)
    D <- cbind(rep(1,n),D)
    A <- diag(RF_weight[i,])
    loc_weight[i,] <- (ginv(t(D)%*%A%*%D)%*%t(D)%*%A)[1,]
    if (Y$type == "distribution"){loc_weight[i,] <- loc_weight[i,]/sum(loc_weight[i,])}
  }
  if (Y$type == "distribution"){
    Res5 <- loc_weight %*% Y_obs
    res5 <- matrix(0,n_xout,nqSup)
    m <- nqSup
    A <- cbind(diag(m), rep(0,m)) + cbind(rep(0,m), -diag(m))
    A <- A[,-c(1,ncol(A))]
    b0 <- rep(0,m-1)
    Pmat <- as(diag(m), "sparseMatrix")
    Amat <- as(t(A), "sparseMatrix")
    for (i in 1:n_xout) {
      temp <- do.call(osqp::solve_osqp,
                    list(P=Pmat, q= -Res5[i,], A=Amat, l=b0, pars = osqp::osqpSettings(verbose = FALSE)))
      res5[i,] <- sort(temp$x)
    }}
  if (Y$type == "covariance"){
    res5 <- array(0,c(m,m,n_xout))
    for (i in 1:n_xout) {
      res5[,,i] <- estcov(M_obs , method=method,weights=loc_weight[i,] ,alpha=alpha)$mean
    }}
  fail <- 0
  if (Y$type == "sphere"){
    #converting 3d coordinates to latitude and longitude coordinates to calculate weighted mean
    S2_obs <- matrix(0,n,2)
    for (i in 1:n){
      S2_obs[i,] = Trans.sph(S_obs[i,]) 
    }
    Res5 <- matrix(0,n_xout,2)
    for (i in 1:n_xout) {
      Res5[i,] <- IntrinsicMean(S2_obs, weights = loc_weight[i,])
      if (Res5[i,1]==5) {fail <- 1;break}
    }
    if (fail==0) {
      #convert back to 3d coordinates
      res5 <- matrix(0,n_xout,3)
      for (i in 1:n_xout) {
        res5[i,] <- Trans.Euclid(Res5[i,])
      }
    }}
  if (fail==0) {return(list(res=res5, weight=loc_weight))}
  if (fail==1) {return(1)}
}


#######################matrix exp############################
matrix.exp <- function(A){
  eig <- eigen(A)
  EA=eig$vectors%*%diag(exp(eig$values))%*%t(eig$vectors)
  return((EA+t(EA))/2)
}

#######################matrix log############################
matrix.log <- function(A){
  eig <- eigen(A)
  LA=eig$vectors%*%diag(log(eig$values))%*%t(eig$vectors)
  return((LA+t(LA))/2)
}


