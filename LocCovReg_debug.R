LocCovReg= function(x,y=NULL,M=NULL,xout,optns = list()){
  if(is.null(optns$metric)){
    metric <- "frobenius"
  } else {
    metric <- optns$metric
  }
  if(!metric %in% c("frobenius","power","cholesky","log_cholesky")){
    stop("metric choice not supported.")
  }
  if(metric=="frobenius") {
    res <- LFRCov(x=x, y=y,M=M,xout=xout,optns = optns)
  } else if(metric=="power") {
    res <- LFRCovPower(x=x, y=y,M=M,xout=xout,optns = optns)
  } else{
    if (is.null(M))
      stop("M must be input for Cholesky and log-Cholesky metrics; y does not apply.")
    res <- LFRCovCholesky(x=x, M=M, xout=xout, optns = optns)
  }
  class(res) <- "covReg"
  return(res)
}

LFRCov  = function(x, y=NULL,M=NULL, xout,optns = list()){
  if(is.null(optns$corrOut)){
    corrOut=FALSE
  } else{
    corrOut=optns$corrOut
  }
  if(is.null(optns$kernel)){
    kernel = 'gauss'
  } else{
    kernel=optns$kernel
  }
  if(is.null(optns$bwMean)){
    bwMean = NA
  } else{
    bwMean=optns$bwMean
  }
  if(is.null(optns$bwCov)){
    bwCov=NA
  } else{
    bwCov=optns$bwCov
  }
  bw2=bwCov
  
  if(!is.matrix(x)&!is.vector(x)){
    stop('x must be a matrix or vector')
  }
  if(is.vector(x)){
    x<- matrix(x,length(x))
  }
  if(is.vector(xout)){
    xout<- matrix(xout,length(xout))
  }
  
  if(!is.matrix(x)){
    stop('x must be a matrix')
  }
  if(!is.matrix(xout)){
    stop('xout must be a matrix')
  }
  if(ncol(x) != ncol(xout)){
    stop('x and xout must have the same number of columns')
  }
  p = ncol(x)
  if(p > 2){
    stop("The number of dimensions of the predictor x is greater than 2.")
  }
  m = nrow(xout)
  
  Kern=kerFctn(kernel)
  K = function(x,h){
    k = 1
    for(i in 1:p){
      k=k*Kern(x[,i]/h[i])
    }
    return(as.numeric(k))
  }
  
  computeLFR=function(idx,x0,bw2){
    #x0 and bw2 are in R^p
    x=as.matrix(x[idx,])
    aux=K(x-matrix(t(x0),nrow=length(idx),ncol=length(x0),byrow=TRUE),bw2) #/ prod(bw2)
    mu0 = mean(aux)
    mu1 = colMeans(aux*(x - matrix(t(x0),nrow=length(idx),ncol=length(x0),byrow=TRUE)))
    mu2=0
    for(i in 1:length(idx)){
      mu2 = mu2 + aux[i]*(x[i,]-x0) %*% t(x[i,]-x0)/length(idx)
    }
    sL = array(0,length(idx))
    for(i in 1:length(idx)){
      sL[i] =aux[i]*(1-t(mu1)%*%solve(mu2)%*%(x[i,]-x0))
    }
    s = sum(sL)
    if(s == 0){
      stop('Bandwidth too small1')
    }
    
    M_aux=array(0,c(dim(M)[1],dim(M)[1],1))
    for(i in 1:length(idx)){
      M_aux[,,1]=M_aux[,,1]+sL[i]*M[,,idx[i]]/s
    }
    M_aux[,,1]
  }
  
  if(!is.null(y)){
    if(!is.matrix(y)){
      stop('y must be a matrix')
    }
    if(nrow(x) != nrow(y)){
      stop('x and y must have the same number of rows')
    }
    n = nrow(y)
    nGrid = ncol(y)
    cm = mean4LocCovReg(x=x,y=y,xout=x,list(bwMean = bwMean))
    bwMean = cm$optns$bwMean
    cmh = cm$mean_out
    
    M=array(0,c(dim(y)[2], dim(y)[2], dim(y)[1]))
    for(i in 1:n){
      M[,,i] = (y[i,] - cmh[i,]) %*% t(y[i,] - cmh[i,])
    }
    if(is.na(sum(bw2))){
      bw2 = bwMean
    }
  } else{
    if(!is.null(M)){
      if(is.list(M)){
        M=array(as.numeric(unlist(M)), dim=c(dim(M[[1]])[1],dim(M[[1]])[1],length(M)))
      }else{
        if(!is.array(M)){
          stop('M must be an array or a list')
        }
      }
      if(nrow(x)!=dim(M)[3]){
        stop("The number of rows of x must be the same as the number of covariance matrices in M")
      }
      #CV for bw2 selection
      if(is.na(sum(bw2))){
        if(p==1){
          bw_choice=SetBwRange(as.vector(x), as.vector(xout), kernel)
          objF=matrix(0,nrow=20,ncol=1)
          aux1=as.matrix(seq(bw_choice$min,bw_choice$max,length.out=20))
          for(i in 1:20){
            for(j in 1:dim(x)[1]){
              aux=as.matrix(Matrix::nearPD(computeLFR(setdiff(1:dim(x)[1],j),x[j],aux1[i]),corr = FALSE)$mat)-M[,,j]
              objF[i]=objF[i]+sum(diag(aux%*%t(aux)))
            }
          }
          ind=which(objF==min(objF))[1]
          bwCV=aux1[ind]
        }
        if(p==2){
          bw_choice1=SetBwRange(as.vector(x[,1]), as.vector(xout[,1]), kernel)
          bw_choice2=SetBwRange(as.vector(x[,2]), as.vector(xout[,2]), kernel)
          objF=matrix(0,nrow=10,ncol=10)
          aux1=seq(bw_choice1$min,bw_choice1$max,length.out=10)
          aux2=seq(bw_choice2$min,bw_choice2$max,length.out=10)
          for(i1 in 1:10){
            for(i2 in 1:10){
              for(j in 1:dim(x)[1]){
                aux=as.matrix(Matrix::nearPD(computeLFR(setdiff(1:dim(x)[1],j),x[j,],c(aux1[i1],aux2[i2])),corr = FALSE)$mat)-M[,,j]
                objF[i1,i2]=objF[i1,i2]+sum(diag(aux%*%t(aux)))
              }
            }
          }
          ind=which(objF==min(objF),arr.ind = TRUE)
          bwCV=c(aux1[ind[1]],aux2[ind[2]])
        }
        bw2=bwCV
      }
    } else{
      stop("y or M must be provided.")
    }
  }
  
  Mout = list()
  if(corrOut){
    for(j in 1:m){
      x0 = xout[j,]
      aux=as.matrix(Matrix::nearPD(computeLFR(1:dim(x)[1],x0,bw2),corr=FALSE)$mat)
      D=diag(1/sqrt(diag(aux)))
      aux=D%*%aux%*%D
      aux=as.matrix(Matrix::forceSymmetric(aux))
      Mout = c(Mout,list(aux))
    }
  } else{
    for(j in 1:m){
      x0 = xout[j,]
      Mout = c(Mout,list(as.matrix(Matrix::nearPD(computeLFR(1:dim(x)[1],x0,bw2),corr = FALSE)$mat)))
    }
  }
  optns$corrOut=corrOut
  optns$kernel=kernel
  optns$bwMean=bwMean
  optns$bwCov=bw2
  return(list(xout=xout, Mout=Mout, optns=optns))
}

LFRCovPower= function(x,y=NULL,M=NULL, xout,optns = list()){
  if(is.null(optns$corrOut)){
    corrOut=FALSE
  } else{
    corrOut=optns$corrOut
  }
  if(is.null(optns$kernel)){
    kernel = 'gauss'
  } else{
    kernel=optns$kernel
  }
  if(is.null(optns$bwMean)){
    bwMean = NA
  } else{
    bwMean=optns$bwMean
  }
  bw1=bwMean
  
  if(is.null(optns$bwCov)){
    bwCov=NA
  } else{
    bwCov=optns$bwCov
  }
  bw=bwCov
  
  if(is.null(optns$alpha)){
    alpha=1
  } else{
    alpha=optns$alpha
  }
  if(alpha<0){
    stop('alpha must be non-negative')
  }
  
  if(!is.matrix(x)&!is.vector(x)){
    stop('x must be a matrix or vector')
  }
  if(is.vector(x)){
    x<- matrix(x,length(x))
  }
  if(is.vector(xout)){
    xout<- matrix(xout,length(xout))
  }
  
  if(!is.matrix(x)){
    stop('x must be a matrix')
  }
  if(!is.matrix(xout)){
    stop('xout must be a matrix')
  }
  if(ncol(x) != ncol(xout)){
    stop('x and xout must have the same number of columns')
  }
  if(!is.na(sum(bw))){
    if(sum(bw<=0)>0){
      stop("bandwidth must be positive")
    }
  }
  p = ncol(x)
  if(p>2){
    stop("The number of dimensions of the Euclidean predictor x must be at most 2")
  }
  m = nrow(xout)
  
  Kern=kerFctn(kernel)
  K = function(x,h){
    k = 1
    for(i in 1:p){
      k=k*Kern(x[,i]/h[i])
    }
    return(as.numeric(k))
  }
  
  computeLFR_originalSpace=function(idx,x0,bw2){
    #x0 and bw2 are in R^p
    x=as.matrix(x[idx,])
    aux=K(x-matrix(t(x0),nrow=length(idx),ncol=length(x0),byrow=TRUE),bw2)
    mu0 = mean(aux)
    mu1 = colMeans(aux*(x - matrix(t(x0),nrow=length(idx),ncol=length(x0),byrow=TRUE)))
    mu2=0
    for(i in 1:length(idx)){
      mu2 = mu2 + aux[i]*(x[i,]-x0) %*% t(x[i,]-x0)/length(idx)
    }
    sL = array(0,length(idx))
    for(i in 1:length(idx)){
      sL[i] =aux[i]*(1-t(mu1)%*%solve(mu2)%*%(x[i,]-x0))
    }
    s = sum(sL)
    if(s == 0){
      stop('Bandwidth is too small')
    }
    
    M_hat=array(0,c(dim(M)[1],dim(M)[1],1))
    if(alpha>0){
      for(i in 1:length(idx)){
        P=eigen(M[,,idx[i]])$vectors
        Lambd_alpha=diag(pmax(0,eigen(M[,,idx[i]])$values)**alpha)
        M_alpha=P%*%Lambd_alpha%*%t(P)
        M_hat[,,1]=M_hat[,,1]+sL[i]*M_alpha/s
      }
      M_hat[,,1]=as.matrix(Matrix::nearPD(M_hat[,,1],corr = FALSE)$mat)
      P=eigen(M_hat[,,1])$vectors
      Lambd_alpha=diag(pmax(0,eigen(M_hat[,,1])$values)**(1/alpha))
      M_hat[,,1]=P%*%Lambd_alpha%*%t(P)
      M_hat[,,1]=as.matrix(Matrix::forceSymmetric(M_hat[,,1]))
    } else{
      for(i in 1:length(idx)){
        P=eigen(M[,,idx[i]])$vectors
        Lambd_alpha=diag(log(pmax(1e-30,eigen(M[,,idx[i]])$values)))
        M_alpha=P%*%Lambd_alpha%*%t(P)
        M_hat[,,1]=M_hat[,,1]+sL[i]*M_alpha/s
      }
      M_hat[,,1]=as.matrix(Matrix::nearPD(M_hat[,,1],corr = FALSE)$mat)
      P=eigen(M_hat[,,1])$vectors
      Lambd_alpha=diag(exp(pmax(0,eigen(M_hat[,,1])$values)))
      M_hat[,,1]=P%*%Lambd_alpha%*%t(P)
      M_hat[,,1]=as.matrix(Matrix::forceSymmetric(M_hat[,,1]))
    }
    M_hat[,,1]
  }
  
  if(!is.null(y)){
    if(!is.matrix(y)){
      stop('y must be a matrix')
    }
    if(nrow(x) != nrow(y)){
      stop('x and y must have the same number of rows')
    }
    n = nrow(y)
    nGrid = ncol(y)
    cm = mean4LocCovReg(x=x,y=y,xout=x,optns=list(bwMean = bw1))
    bw1 = cm$optns$bwMean
    cmh = cm$mean_out
    
    M=array(0,c(dim(y)[2], dim(y)[2], dim(y)[1]))
    for(i in 1:n){
      M[,,i] = (y[i,] - cmh[i,]) %*% t(y[i,] - cmh[i,])
    }
    if(is.na(sum(bw))){
      bw = bw1
    }
  } else{
    if(is.null(M)){
      stop("y or M must be provided")
    }
    if(is.list(M)){
      M=array(as.numeric(unlist(M)), dim=c(dim(M[[1]])[1],dim(M[[1]])[1],length(M)))
    } else{
      if(!is.array(M)){
        stop('M must be an array or a list')
      }
    }
    if(nrow(x)!=dim(M)[3]){
      stop("The number of rows of x must be the same as the number of covariance matrices in M")
    }
    
    #CV for bw selection
    if(is.na(sum(bw))){
      if(p==1){
        bw_choice=SetBwRange(as.vector(x), as.vector(xout), kernel)
        objF=matrix(0,nrow=20,ncol=1)
        aux1=as.matrix(seq(bw_choice$min,bw_choice$max,length.out=20))
        for(i in 1:20){
          for(j in 1:dim(x)[1]){
            aux=computeLFR_originalSpace(setdiff(1:dim(x)[1],j),x[j],aux1[i])-M[,,j]
            objF[i]=objF[i]+sum(diag(aux%*%t(aux)))
          }
        }
        ind=which(objF==min(objF))[1]
        bwCV=aux1[ind]
      }
      if(p==2){
        bw_choice1=SetBwRange(as.vector(x[,1]), as.vector(xout[,1]), kernel)
        bw_choice2=SetBwRange(as.vector(x[,2]), as.vector(xout[,2]), kernel)
        if(n<=30){
          objF=matrix(0,nrow=6,ncol=6)
          aux1=seq(bw_choice1$min,bw_choice1$max,length.out=6)
          aux2=seq(bw_choice2$min,bw_choice2$max,length.out=6)
          for(i1 in 1:6){
            for(i2 in 1:6){
              for(j in 1:dim(x)[1]){
                aux=computeLFR_originalSpace(setdiff(1:dim(x)[1],j),x[j,],c(aux1[i1],aux2[i2]))-M[,,j]
                objF[i1,i2]=objF[i1,i2]+sum(diag(aux%*%t(aux)))
              }
            }
          }
          ind=which(objF==min(objF),arr.ind = TRUE)
          bwCV=c(aux1[ind[1]],aux2[ind[2]])
        } else{
          randIndices=sample(dim(x)[1])
          groupIndices=cut(seq(1,dim(x)[1]),breaks=10,labels=FALSE)
          cv10fold_compute=function(v){
            tmp = 0.001*diag(nrow(M[,,v]))
            state = try({tmp = computeLFR_originalSpace(leaveIn,x[v,],c(aux1[i1],aux2[i2]))}, silent = T)
            aux = tmp-M[,,v]
            sum(diag(aux%*%t(aux)))
          }
          objF=matrix(0,nrow=6,ncol=6)
          aux1=seq(bw_choice1$min,bw_choice1$max,length.out=6)
          aux2=seq(bw_choice2$min,bw_choice2$max,length.out=6)
          for(i1 in 1:6){
            for(i2 in 1:6){
              for(j in 1:10){
                leaveIn=setdiff(1:(dim(x)[1]),randIndices[groupIndices==j])
                objF[i1,i2]=objF[i1,i2]+sum(sapply(randIndices[groupIndices==j],cv10fold_compute))
              }
            }
          }
          ind=which(objF==min(objF),arr.ind = TRUE)
          bwCV=c(aux1[ind[1]],aux2[ind[2]])
        }
      }
      bw=bwCV
    }
  }
  Mout = list()
  if(corrOut){
    for(j in 1:m){
      x0 = xout[j,]
      aux=computeLFR_originalSpace(1:dim(x)[1],x0,bw)
      D=diag(1/sqrt(diag(aux)))
      aux=D%*%aux%*%D
      aux=as.matrix(Matrix::forceSymmetric(aux))
      Mout = c(Mout,list(aux))
    }
  } else{
    for(j in 1:m){
      x0 = xout[j,]
      Mout = c(Mout,list(computeLFR_originalSpace(1:dim(x)[1],x0,bw)))
    }
  }
  optns$corrOut=corrOut
  optns$kernel=kernel
  optns$bwMean=bw1
  optns$bwCov=bw
  return(list(xout=xout, Mout=Mout, optns=optns))
}

LFRCovCholesky <- function(x, M, xout, optns=list()){
  if(!is.matrix(x)&!is.vector(x)){
    stop('x must be a matrix or vector')
  }
  if(!is.matrix(xout)&!is.vector(xout)){
    stop('xout must be a matrix or vector')
  }
  
  if(is.vector(x)){x<- matrix(x,length(x)) }
  if(is.vector(xout)){xout<- matrix(xout,length(xout))}
  
  if(ncol(x) != ncol(xout)){
    stop('x and xout must have same number of columns')
  }
  
  
  if(is.null(optns$bwCov)){
    bwCov = NA
  } else {
    bwCov = optns$bwCov
    if(min(bwCov)<=0){
      stop("bandwidth must be positive")
    }
  }
  
  if(is.null(optns$kernel)){
    kernel= 'gauss'
  } else {
    kernel = optns$kernel
  } 
  
  if(is.null(optns$corrOut)){
    corrOut = FALSE
  } else {
    corrOut = optns$corrOut
  }
  
  if(is.null(optns$metric)){
    metric = 'log_cholesky'
  } else {
    metric =  optns$metric
  }
  
  p = ncol(x)
  if(p>2){
    stop("The number of dimensions of the predictor x must be at most 2")
  }
  m = nrow(xout)
  n = nrow(x)
  
  
  Kern=kerFctn(kernel)
  K = function(x,h){
    k = 1
    for(i in 1:p){
      k=k*Kern(x[,i]/h[i])
    }
    return(as.numeric(k))
  }
  if(is.null(M)){
    stop("M must be provided")
  }
  
  if(class(M) == 'array'){
    MM = list()
    if(class(M) == 'array'){
      for (i in 1:dim(M)[3]) {
        MM[[i]] = M[,,i]
      }
    }
    M = lapply(MM, function(X) (X+t(X))/2)
  }else{
    if(!class(M)=="list"){
      stop('M must be an array or a list')
    }
    M = lapply(M, function(X) (X+t(X))/2)
  }
  
  if(nrow(x)!= length(M)){
    stop("the number of rows of x must be the same as the number of covariance matrices in M")
  }
  
  computeLFRSPD=function(idx,x0,bw2){
    #idx: index for x
    #x0 m-by-p matrix,
    #bw2 are in b-by-p
    x=as.matrix(x[idx,])
    aux=K(x-matrix(t(x0),nrow=length(idx),ncol=length(x0),byrow=TRUE),bw2) #/ prod(bw2)
    mu0 = mean(aux)
    mu1 = colMeans(aux*(x - matrix(t(x0),nrow=length(idx),ncol=length(x0),byrow=TRUE)))
    mu2=0
    for(i in 1:length(idx)){
      mu2 = mu2 + aux[i]*(x[i,]-x0) %*% t(x[i,]-x0)/length(idx)
    }
    sL = array(0,length(idx))
    # browser()
    for(i in 1:length(idx)){
      # state = try(solve(mu2))
      # if(class(state) == "try-error") browser()
      tmp = 0
      try({tmp = aux[i]*(1-t(mu1)%*%solve(mu2)%*%(x[i,]-x0))}, silent = T)
      sL[i] = tmp
    }
    s = sum(sL)
    if(s == 0){
      stop('Bandwidth too small2')
    }
    
    Mout = list()
    MM = M[idx]
    n = length(idx)
    if(metric == 'log_cholesky'){
      LL = lapply(MM, chol)
      L = lapply(LL, function(X) X - diag(diag(X)))
      D = lapply(LL, function(X) diag(X))
      
      U = 0
      E = 0
      for (i in 1:n) {
        U = U + sL[i]*L[[i]]
        E = E + sL[i]*log(D[[i]])
      }
      SS = U/s + diag(exp(E/s))
      Mout = t(SS)%*%SS
      
    } else {
      L = lapply(MM, chol)
      U = 0
      for (i in 1:n) {
        U = U + sL[i]*L[[i]]
      }
      Mout = t(U/s) %*% (U/s)
    }
    
    return(Mout)
  }
  
  distance <- function(M1, M2){
    if(metric == 'log_cholesky'){
      LL1 = chol(M1); LL2 = chol(M2)
      L1 = LL1 - diag(diag(LL1)); L2 = LL2 - diag(diag(LL2))
      D1 = diag(LL1); D2 = diag(LL2)
      L = L1 - L2; D = log(D1) - log(D2)
      res = sqrt(sum(sum(L^2))+sum(D^2))
    }else{ 
      L1 = chol(M1); L2 = chol(M2)
      L = L1 - L2;
      res = sqrt(sum(sum(L^2)))
    }
    return(res)
  }
  
  #CV for bwCov selection
  if(is.na(sum(bwCov))){
    if(p==1){
      bw_choice=SetBwRange(as.vector(x), as.vector(xout), kernel)
      objF=matrix(0,nrow=20,ncol=1)
      aux1=as.matrix(seq(bw_choice$min,bw_choice$max,length.out=20))
      for(i in 1:20){
        for(j in 1:dim(x)[1]){
          distj = distance(computeLFRSPD(setdiff(1:dim(x)[1],j),x[j],aux1[i]), M[[j]])
          objF[i]=objF[i] + distj
        }
      }
      ind=which(objF==min(objF))[1]
      bwCV=aux1[ind]
    }
    if(p==2){
      bw_choice1=SetBwRange(as.vector(x[,1]), as.vector(xout[,1]), kernel)
      bw_choice2=SetBwRange(as.vector(x[,2]), as.vector(xout[,2]), kernel)
      if(n<=30){
        objF=matrix(0,nrow=6,ncol=6)
        aux1=seq(bw_choice1$min,bw_choice1$max,length.out=6)
        aux2=seq(bw_choice2$min,bw_choice2$max,length.out=6)
        for(i1 in 1:6){
          for(i2 in 1:6){
            for(j in 1:dim(x)[1]){
              distj=distance(computeLFRSPD(setdiff(1:dim(x)[1],j),x[j,],c(aux1[i1],aux2[i2])), M[[j]])
              objF[i1,i2]=objF[i1,i2]+distj
            }
          }
        }
        ind=which(objF==min(objF),arr.ind = TRUE)
        bwCV=c(aux1[ind[1]],aux2[ind[2]])
      }else{
        randIndices=sample(dim(x)[1])
        groupIndices=cut(seq(1,dim(x)[1]),breaks=10,labels=FALSE)
        cv10fold_compute=function(v){
          tmp = 0.001*diag(nrow(M[[v]]))
          state = try({tmp = computeLFRSPD(leaveIn,x[v,],c(aux1[i1],aux2[i2]))}, silent = T)
          # if(class(state) == "try-error") browser()
          distance(tmp,M[[v]])
        }
        objF=matrix(0,nrow=6,ncol=6)
        aux1=seq(bw_choice1$min,bw_choice1$max,length.out=6)
        aux2=seq(bw_choice2$min,bw_choice2$max,length.out=6)
        for(i1 in 1:6){
          for(i2 in 1:6){
            for(j in 1:10){
              leaveIn=setdiff(1:(dim(x)[1]),randIndices[groupIndices==j])
              objF[i1,i2]=objF[i1,i2]+sum(sapply(randIndices[groupIndices==j],cv10fold_compute))
            }
          }
        }
        ind=which(objF==min(objF),arr.ind = TRUE)
        bwCV=c(aux1[ind[1]],aux2[ind[2]])
      }
    }
    bwCov=bwCV
  }
  
  
  Mout = list()
  for (j in 1:nrow(xout)) {
    Mout[[j]] = computeLFRSPD(1:dim(x)[1], xout[j,], bwCov)
  }
  
  if(corrOut){
    for(j in 1:nrow(xout)){
      D=diag(1/sqrt(diag(Mout[[j]])))
      Mout[[j]]=D%*%Mout[[j]]%*%D
      Mout[[j]]=as.matrix(Matrix::forceSymmetric(Mout[[j]]))
    }
  }
  out = list(xout=xout,Mout=Mout,optns=list(bwCov =bwCov,kernel=kernel,corrOut=corrOut,metric=metric))
  return(out)
}

kerFctn <- function(kernel_type){
  if (kernel_type=='gauss'){
    ker <- function(x){
      dnorm(x) #exp(-x^2 / 2) / sqrt(2*pi)
    }
  } else if(kernel_type=='rect'){
    ker <- function(x){
      as.numeric((x<=1) & (x>=-1))
    }
  } else if(kernel_type=='epan'){
    ker <- function(x){
      n <- 1
      (2*n+1) / (4*n) * (1-x^(2*n)) * (abs(x)<=1)
    }
  } else if(kernel_type=='gausvar'){
    ker <- function(x) {
      dnorm(x)*(1.25-0.25*x^2)
    }
  } else if(kernel_type=='quar'){
    ker <- function(x) {
      (15/16)*(1-x^2)^2 * (abs(x)<=1)
    }
  } else {
    stop('Unavailable kernel')
  }
  return(ker)
}

SetBwRange <- function(xin, xout, kernel_type) {
  xinSt <- sort(xin)
  bw.min <- max(diff(xinSt), xinSt[2] - min(xout), max(xout) - xinSt[length(xin)-1])*1.1 / (ifelse(kernel_type == "gauss", 3, 1) * ifelse(kernel_type == "gausvar", 2.5, 1))
  bw.max <- diff(range(xin))/3 / (ifelse(kernel_type == "gauss", 3, 1) * ifelse(kernel_type == "gausvar", 2.5, 1))
  if (bw.max < bw.min) {
    if (bw.min > bw.max*3/2) {
      #warning("Data is too sparse.")
      bw.max <- bw.min*1.01
    } else bw.max <- bw.max*3/2
  }
  return(list(min=bw.min, max = bw.max))
}

bwCV <- function(xin, qin, xout, optns) {
  p=ncol(xin)
  if(p==1){
    compareRange <- (xin[,1] > min(xin[,1]) + diff(range(xin[,1]))/5) & (xin[,1] < max(xin[,1]) - diff(range(xin[,1]))/5)
  }else{
    compareRange <- (xin[,1] > min(xin[,1]) + diff(range(xin[,1]))/5) & (xin[,1] < max(xin[,1]) - diff(range(xin[,1]))/5) & (xin[,2] > min(xin[,2]) + diff(range(xin[,2]))/5) & (xin[,2] < max(xin[,2]) - diff(range(xin[,2]))/5)
  }
  
  # k-fold
  objFctn <- function(bw) {
    optns1 <- optns
    optns1$bw <- bw
    folds <- numeric(nrow(xin))
    nn <- sum(compareRange)
    numFolds <- ifelse(nn > 30, 10, sum(compareRange))
    
    tmp <- c(sapply(1:ceiling(nn/numFolds), function(i)
      sample(x = seq_len(numFolds), size = numFolds, replace = FALSE)))
    tmp <- tmp[1:nn]
    repIdx <- which(diff(tmp) == 0)
    for (i in which(diff(tmp) == 0)) {
      s <- tmp[i]
      tmp[i] <- tmp[i-1]
      tmp[i-1] <- s
    }
    #tmp <- cut(1:n,breaks = seq(0,n,length.out = numFolds+1), labels=FALSE)
    #tmp <- tmp[sample(seq_len(n), n)]
    
    folds[compareRange] <- tmp
    
    qfit <- lapply(seq_len(numFolds), function(foldidx) {
      testidx <- which(folds == foldidx)
      res <- LocWassReg(xin = xin[-testidx,], qin = qin[-testidx,], xout = xin[testidx,],
                        optns = optns1)
      res # each row is a qt function
    })
    qfit <- do.call(rbind, qfit)
    mean(apply((qfit - qin[which(compareRange)[order(tmp)],])^2, 1, pracma::trapz, x = optns1$qSup))
  }
  
  if(p==1){
    aux=SetBwRange(xin = xin[,1], xout = xout[,1], kernel_type = optns$ker)
    bwRange <- matrix(c(aux$min,aux$max),nrow=2,ncol=1)
  }else{
    aux=SetBwRange(xin = xin[,1], xout = xout[,1], kernel_type = optns$ker)
    aux2=SetBwRange(xin = xin[,2], xout = xout[,2], kernel_type = optns$ker)
    bwRange <- as.matrix(cbind(c(aux$min,aux$max),c(aux2$min,aux2$max)))
  }
  if(!is.null(optns$bwRange)){
    if(p==1){
      if (min(optns$bwRange) < min(bwRange)) {
        message("Minimum bandwidth is too small and has been reset.")
      }else{
        bwRange[1,1] <- min(optns$bwRange)
      }
      if (max(optns$bwRange) >  min(bwRange)) {
        bwRange[2,1] <- max(optns$bwRange)
      }else {
        message("Maximum bandwidth is too small and has been reset.")
      }
    }else{
      #Check for first dimension of the predictor
      if (min(optns$bwRange[,1]) < min(bwRange[,1])) {
        message("Minimum bandwidth of first predictor dimension is too small and has been reset.")
      }else{
        bwRange[1,1] <- min(optns$bwRange[,1])
      }
      if (max(optns$bwRange[,1]) >  min(bwRange[,1])) {
        bwRange[2,1] <- max(optns$bwRange[,1])
      } else {
        message("Maximum bandwidth of first predictor dimension is too small and has been reset.")
      }
      #Check for second dimension of the predictor
      if (min(optns$bwRange[,2]) < min(bwRange[,2])) {
        message("Minimum bandwidth of second predictor dimension is too small and has been reset.")
      }else{
        bwRange[1,2] <- min(optns$bwRange[,2])
      }
      if (max(optns$bwRange[,2]) >  min(bwRange[,2])) {
        bwRange[2,2] <- max(optns$bwRange[,2])
      }else{
        message("Maximum bandwidth of second predictor dimension is too small and has been reset.")
      }
    }
  }
  if(p==1){
    res <- optimize(f = objFctn, interval = bwRange[,1])$minimum
  }else{
    res <- optim(par=rowMeans(bwRange),fn=objFctn,lower=bwRange[,1],upper=bwRange[,2],method='L-BFGS-B')$par
  }
  res
}

