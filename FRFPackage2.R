##Paper:Random Forest Weighted Local Fr\'echet Regression with Random Objects#########
##This R package is developed from the FrechForest package, with just a few more regression scenarios (distributions, spd and spherical data) added to the content. We sincerely thank the original author Louis Capitaine. See “https://github.com/Lcapitaine/FrechForest” for FrechForest package.




#Frechet random forest functions

#' Read the parameters of the function
#'
#' @param z
#'
#'
#' @keywords internal
read.Xarg <- function(z){
  type <- NULL
  issou <- rep(NA, length(z))
  for (i in 1:length(z)){
    issou[i] <- is.null(z[i])
    if (issou[i]==FALSE) type <- c(type,z[i]$type)
  }
  return(type)
}

#' Factor partitions finder
#'
#' This function is used to find all the unique partitions of k factors into 2 groups
#'
#' @param Factor
#' @param id
#'
#' @keywords internal
Fact.partitions <- function(Factor, id){
  
  U <- unique(Factor)
  P <- Part.facts[[length(U)]]
  L <- list()
  for (k in 1:nrow(P)){
    w <- which(P[k,]==0)
    U_courant <- U[w]
    W <- NULL
    for (m in U_courant){
      W <- c(W,which(Factor==m))
    }
    L[[k]] <- id[W]
  }
  return(L)
}

#' Ordonne
#'
#' @param X
#' @param time
#' @param id
#'
#'
#' @keywords internal
ordonne <- function(X , time , id){
  mat  <- matrix(NA, length(unique(id)), length(unique(time)))
  for( i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    t_w <- time[w]
    w_time <- NULL
    for (j in 1:length(w)){
      w_time <- c(w_time, which(unique(time)==t_w[j]))
    }
    mat[i,w_time] <- X[w]
  }
  return(mat)
}

#' Impurity
#'
#' Compute the impurity of a given vector
#'
#' @param Y
#' @param timeScale
#'
#' @import kmlShape
#' @import Evomorph
#' @import geomorph
#' @import stats
#'
#'
#' @keywords internal
impurity <- function(Y, timeScale=0.1){
  
  if (Y$type=="curve"){
    traj <- Y$Y
    id <- Y$id
    time <- Y$time
    imp <- 0
    trajLong <- data.frame(id=id,time=time,traj=traj)
    meanF <- meanFrechet(trajLong = trajLong, timeScale = timeScale)
    for (i in unique(id)){
      imp <- imp + distFrechet(meanF$times, meanF$traj, time[which(id==i)], traj[which(id==i)], timeScale = timeScale)^2
    }
    imp <- imp/length(unique(id))
    return(imp)
  }
  
  if (Y$type=="image"){
    if (length(Y$id)==1){
      return(0)
    }
    return(mean(apply(Y$Y,2,"var")))
  }
  
  if (Y$type=="scalar"){
    if (length(Y$Y)==1){
      return(0)
    }
    return(var(Y$Y))
  }
  
  if (Y$type=="factor"){
    p = table(Y$Y)/length(Y$Y)
    return(-sum(p*log2(p)))
  }
  
  if (Y$type=="shape"){
    ms<- mshape(Y$Y[,,,drop=FALSE])
    return(mean(ShapeDist(Y$Y,ms)^2))
  }
  
  if (Y$type=="distribution"){
    if(!is.matrix(Y$Y)) Y$Y = matrix(Y$Y, nrow = 1)
    id <- Y$id
    imp <- 0
    #meanF <- DenFMean(qin = Y$Y, optns = list(qSup = qSup))$qout
    meanF <- colMeans(Y$Y)
    #d1 <- list(x = qSup,y = meanF)
    for (i in unique(id)){
      #d2 <- list(x = qSup,y = Y$Y[which(id==i),])
      #imp <- imp + dist4den(d1 = d1, d2 = d2, fctn_type = "quantile", optns = list(nqSup=nqSup))^2
      imp <- imp + pracma::trapz(qSup, (meanF - Y$Y[which(id==i),])^2)
    }
    imp <- imp/length(unique(id))
    return(imp)
  }
  
  if (Y$type=="sphere"){
    if(!is.matrix(Y$Y)) Y$Y = matrix(Y$Y, nrow = 1)
    id <- Y$id
    imp <- 0
    S <- list()
    nRow <- nrow(Y$Y)
    for (h in 1:nRow){
      S[[h]] = Y$Y[h,] # project onto Sphere
    }
    S = riemfactory(S, name="sphere")  
    meanF <- as.vector(rbase.mean(S)$x)
    for (i in unique(id)){
      tmp = crossprod(meanF,Y$Y[which(id==i),])[1,1]
      tmp = sign(tmp) * min(1, abs(tmp))
      imp <- imp + (acos(tmp))^2
    }
    imp <- imp/length(unique(id))
    return(imp)
  }

  if (Y$type=="sphere2"){
    if(!is.matrix(Y$Y)) Y$Y = matrix(Y$Y, nrow = 1)
    id <- Y$id
    imp <- 0
    meanF <- IntrinsicMean(Y$Y)
    for (i in unique(id)){
      imp <- imp + distCosine(meanF,Y$Y[which(id==i),],r=1)^2
    }
    imp <- imp/length(unique(id))
    return(imp)
  }

  if (Y$type=="covariance"){
    if(is.na(dim(Y$Y)[3])) Y$Y = array(Y$Y, dim = c(dim(Y$Y)[1],dim(Y$Y)[2],1))
    id <- Y$id
    imp <- 0
    meanF <- estcov(Y$Y, method= method,alpha=alpha)$mean
    for (i in unique(id)){
      imp <- imp + distcov(meanF, Y$Y[,,which(id==i)], method=method, alpha=alpha)^2
    }
    imp <- imp/length(unique(id))
    return(imp)
  }

}

#' Impurity Split
#'
#' @param Y
#' @param split
#' @param timeScale
#'
#' @import survival
#' @keywords internal
impurity_split <- function(Y,split,timeScale=0.1){
  impur <- 0
  imp <- list()
  for (i in 1:2){
    fils <- unique(Y$id)[which(split==i)]
    prop <- length(fils)/length(unique(Y$id))
    if (Y$type=="curve"){
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w, which(Y$id==fils[j]))
      }
      imp[[i]] <- impurity(list(type="curve",Y=Y$Y[w],id=Y$id[w],time=Y$time[w]))
      impur <- impur + imp[[i]]*prop
    }
    
    if (Y$type=="image"){
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w,which(Y$id==fils[j]))
      }
      imp[[i]] <- impurity(list(type="image", Y=Y$Y[w,], id=Y$id[w]))
      impur <- impur + imp[[i]]*prop
    }
    
    if (Y$type=="shape"){
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w, which(Y$id==fils[j]))
      }
      if (length(w)>1){imp[[i]] <- impurity(list(type=Y$type,Y=Y$Y[,,w],id=Y$id[w]))
      impur <- impur + imp[[i]]*prop}
      else {imp[[i]] <- 0}
    }
    
    if (Y$type=="scalar" || Y$type=="factor") {
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w, which(Y$id==fils[j]))
      }
      imp[[i]] <- impurity(list(type=Y$type,Y=Y$Y[w],id=Y$id[w]))
      impur <- impur + imp[[i]]*prop
    }
    
    if (Y$type == "surv"){
      
      tryCatch({
        sp=NULL
        for (i in 1:length(split)){
          sp = c(sp, rep(split[i],length(which(Y$id==unique(Y$id)[i]))))
        }
        impur= 1/(1+survdiff(Surv(Y$time,Y$Y)~sp)$chisq)}, error = function(sp){impur=Inf})
      
    }
    
    if (Y$type=="distribution") {
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w, which(Y$id==fils[j]))
      }
      imp[[i]] <- impurity(list(type="distribution",Y=Y$Y[w,],id=Y$id[w]))
      impur <- impur + imp[[i]]*prop
    }

    if (Y$type=="sphere") {
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w, which(Y$id==fils[j]))
      }
      imp[[i]] <- impurity(list(type="sphere",Y=Y$Y[w,],id=Y$id[w]))
      impur <- impur + imp[[i]]*prop
    }

    if (Y$type=="sphere2") {
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w, which(Y$id==fils[j]))
      }
      imp[[i]] <- impurity(list(type="sphere2",Y=Y$Y[w,],id=Y$id[w]))
      impur <- impur + imp[[i]]*prop
    }

    if (Y$type=="covariance") {
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w, which(Y$id==fils[j]))
      }
      imp[[i]] <- impurity(list(type="covariance",Y=Y$Y[,,w],id=Y$id[w]))
      impur <- impur + imp[[i]]*prop
    }
    
  }
  return(list(impur=impur, imp_list=imp))
}

#' Extremely randomized split
#'
#' @param X
#' @param Y
#' @param timeScale
#' @param ntry
#'
#' @import kmlShape
#' @import Evomorph
#'
#' @keywords internal
ERvar_split <- function(X ,Y,ntry=3,timeScale=0.1){
  
  impur <- rep(0,dim(X$X)[length(dim(X$X))])
  toutes_imp <- list()
  split <- list()
  Pure <- FALSE
  
  for (i in 1:dim(X$X)[length(dim(X$X))]){
    
    if (X$type=="factor"){
      
      if (length(unique(X$X[,i]))>1){
        L <- Fact.partitions(X$X[,i],X$id)
        split_courant <- list()
        impur_courant <- rep(NA,length(L))
        toutes_imp_courant <- list()
        
        # On tire une partition au hasard
        tirage <- sample(1:length(L), 1)
        # Il faut maintenant regarder quelles sont les meilleures combinaisons ::
        
        split[[i]] <- rep(2,length(X$id))
        for (l in L[[tirage]]){
          split[[i]][which(X$id==l)] <- 1
        }
        # Il faut maintenant regarder la qualite du decoupage ::
        impurete <- impurity_split(Y,split[[i]])
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }
      else {
        impur[i] <- Inf
        split[[i]] <- Inf
      }
    }
    
    if( X$type=="curve"){
      
      # Il faut commencer par tirer les multiples centres ::
      
      id_centers <- matrix(NA,ntry,2)
      for (l in 1:ntry){
        id_centers[l,] <- sample(unique(X$id),2)
      }
      
      ### Il faut ensuite boucler sur le ntry
      split_prime <- matrix(2,ntry,length(unique(X$id)))
      u <- 0
      impurete2 <- list()
      qui <- NULL
      imp <- NULL
      
      for (c in 1:ntry){
        
        w_gauche <- which(X$id==id_centers[c,1])
        w_droit <- which(X$id==id_centers[c,2])
        
        for (l in 1:length(unique(X$id))){
          
          w <- which(X$id==unique(X$id)[l])
          dg <- distFrechet(X$time[w_gauche],X$X[w_gauche,i],X$time[w],X$X[w,i], timeScale = timeScale)
          dd <- distFrechet(X$time[w_droit],X$X[w_droit,i],X$time[w],X$X[w,i], timeScale = timeScale)
          if (dg<=dd) split_prime[c,l] <- 1
        }
        
        if (length(unique(split_prime[c,]))>1){
          u <- u+1
          qui <- c(qui, c)
          impurete2[[c]] <- impurity_split(Y,split_prime[c,], timeScale)
          imp <- c(imp,impurete2[[c]]$impur)
        }
        
      }
      
      if (u>0){
        gagnant <- qui[which.min(imp)]
        split[[i]] <- split_prime[gagnant,]
        impurete <- impurete2[[gagnant]]
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }
      
      else{
        impur[i] <- Inf
        split[[i]] <- Inf
      }
      
    }
    
    if (X$type=="shape"){
      if (dim(X$X[,,,i])[3]>2){
        
        id_centers <- matrix(NA,ntry,2)
        for (l in 1:ntry){
          id_centers[l,] <- sample(X$id,2)
        }
        
        split_prime <- matrix(2,ntry,length(X$id))
        
        u <- 0
        qui <- NULL
        impurete2 <- list()
        imp <- NULL
        
        for (c in 1:ntry){
          dg <- ShapeDist(X$X[,,,i],X$X[,,which(X$id==id_centers[c,1]),i])
          dd <- ShapeDist(X$X[,,,i],X$X[,,which(X$id==id_centers[c,2]),i])
          for (l in 1:length(unique(X$id))){
            if (dg[l]<=dd[l]) split_prime[c,l] <- 1
          }
          if (length(split_prime[c,])>1){
            u <- u+1
            qui <- c(qui,c)
            impurete2[[c]] <- impurity_split(Y,split_prime[c,], timeScale)
            imp <- c(imp, impurete2[[c]]$impur)
          }
        }
        ## Il nous faut calculer les distances à gauche et à droite pour chaque element
        
        if (u>0){
          gagnant <- qui[which.min(imp)]
          split[[i]] <- split_prime[gagnant,]
          
          impurete <- impurete2[[gagnant]]
          impur[i] <- impurete$impur
          toutes_imp[[i]] <- impurete$imp_list
        }
        
        else{
          impur[i] <- Inf
          split[[i]] <- Inf
        }
        
      }
      
      else{
        split[[i]] <- c(1,2)
        impurete <- impurity_split(Y,split[[i]], timeScale)
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }
    }
    
    if (X$type=="image"){
      if (nrow(X$X)>2){
        id_centers <- matrix(NA,ntry,2)
        for (l in 1:ntry){
          id_centers[l,] <- sample(X$id,2)
        }
        
        split_prime <- matrix(2,ntry,length(X$id))
        
        
        u <- 0
        qui <- NULL
        impurete2 <- list()
        imp <- NULL
        
        for (c in 1:ntry){
          
          w_g <- which(X$id==id_centers[c,1])
          w_d <- which(X$id==id_centers[c,2])
          ### Il nous faut calculer la distance :
          dg = apply(apply(X$X[,,i],1,"-",X$X[w_g,,i])^2,2,"mean")
          dd = apply(apply(X$X[,,i],1,"-",X$X[w_d,,i])^2,2,"mean")
          
          split_prime[c,which((dg<=dd)==TRUE)]=1
          if (length(unique(split_prime[c,]))>1){
            u <-u+1
            qui <- c(qui,c)
            impurete2[[c]] <- impurity_split(Y,split_prime[c,], timeScale)
            imp <- c(imp,impurete2[[c]]$impur)
          }
          
        }
        
        
        
        if (u>0){
          gagnant <- qui[which.min(imp)]
          split[[i]] <- split_prime[gagnant,]
          impurete <- impurete2[[gagnant]]
          impur[i] <- impurete$impur
          toutes_imp[[i]] <- impurete$imp_list
        }
        
        else{
          impur[i] <- Inf
          split[[i]] <- Inf
        }
        
      }
      
      else{
        split[[i]] <- c(1,2)
        impurete <- impurity_split(Y,split[[i]], timeScale)
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }
    }
    
    if(X$type=="scalar"){
      if (length(unique(X$X[,i]))>2){
        
        ### On doit tier les centres
        #centers <- sample(X$X[,i],2)
        
        centers <- matrix(NA,ntry,2)
        for (l in 1:ntry){
          centers[l,] <- sample(X$X[,i],2)
        }
        
        #split[[i]] <- rep(2,length(X$X[,i]))
        split_prime <- matrix(2,ntry,length(X$X[,i]))
        
        for (l in 1:length(X$X[,i])){
          for (k in 1:ntry){
            if (abs(centers[k,1]-X$X[l,i])<= abs(centers[k,2]-X$X[l,i])) split_prime[k,l] <- 1
          }
        }
        
        u <- 0
        qui <- NULL
        impurete2 <- list()
        imp <- NULL
        for (k in 1:ntry){
          if (length(unique(split_prime[k,]))>1){
            u <- u+1
            qui <- c(qui,k)
            impurete2[[k]] <- c(impurete2,impurity_split(Y,split_prime[k,], timeScale))
            imp <- c(imp, impurete2[[k]]$impur)
          }
        }
        
        if (u>0){
          gagnant <- qui[which.min(imp)]
          split[[i]] <- split_prime[gagnant,]
          impurete <- impurete2[[gagnant]]
          impur[i] <- impurete$impur
          toutes_imp[[i]] <- impurete$imp_list
        }
        
        else{
          impur[i] <- Inf
          split[[i]] <- Inf
        }
      }
      
      else {
        impur[i] <- Inf
        split[[i]] <- Inf
      }
    }
  }
  
  if (length(unique(impur))==1 & is.element(Inf,impur)==TRUE){
    return(list(Pure=TRUE))
  }
  true_split <- which.min(impur)
  split <- split[[true_split]]
  return(list(split=split, impurete=min(impur),impur_list = toutes_imp[[true_split]], variable=which.min(impur), Pure=Pure))
}

#' Classical Variable Split function
#'
#' @param X
#' @param Y
#' @param timeScale
#'
#' @import kmlShape
#' @import Evomorph
#'
#' @keywords internal
var_split <- function(X ,Y,timeScale=0.1){
  # Pour le moment on se concentre sur le cas des variables courbes ::
  impur <- rep(0,dim(X$X)[length(dim(X$X))])
  toutes_imp <- list()
  split <- list()
  centers <- list() # On va stocker les centres associes aux kmeans
  Pure <- FALSE
  
  for (i in 1:dim(X$X)[length(dim(X$X))]){
    
    if (X$type=="factor"){
      if (length(unique(X$X[,i]))>1){
        L <- Fact.partitions(X$X[,i],X$id)
        split_courant <- list()
        impur_courant <- rep(NA,length(L))
        toutes_imp_courant <- list()
        # Il faut maintenant regarder quelles sont les meilleures combinaisons ::
        for (k in 1:length(L)){
          split_courant[[k]] <- rep(2,length(X$id))
          for (l in L[[k]]){
            split_courant[[k]][which(X$id==l)] <- 1
          }
          # Il faut maintenant regarder la qualite du decoupage ::
          impurete <- impurity_split(Y,split_courant[[k]])
          impur_courant[k] <- impurete$impur
          toutes_imp_courant[[k]] <- impurete$imp_list
        }
        select <- which.min(impur_courant)
        split[[i]] <- split_courant[[select]]
        impur[i] <- impur_courant[select]
        toutes_imp[[i]] <- toutes_imp_courant[[select]]
      }
      else {
        impur[i] <- Inf
        split[[i]] <- Inf
      }
    }
    
    if( X$type=="curve"){
      mclds <- kmlShape::cldsWide(ordonne(X$X[,i], X$time, X$id), unique(X$time), unique(X$id))
      crit <- kmlShape::kmlShape(mclds, nbClusters = 2, timeScale = timeScale, toPlot="none")
      att <- attributes(crit)
      split[[i]] <- att$clusters
      impurete <- impurity_split(Y,split[[i]], timeScale)
      impur[i] <- impurete$impur
      toutes_imp[[i]] <- impurete$imp_list
    }
    
    if(X$type=="scalar"){
      if (length(unique(X$X[,i]))>2){
        sp <- kmeans(X$X[,i], centers=2)
        split[[i]] <- sp$cluster
        impurete <- impurity_split(Y,split[[i]], timeScale)
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }
      
      if (length(unique(X$X[,i]))==2){
        split[[i]] <- rep(2,length(X$X[,i]))
        split[[i]][which(X$X[,i]==unique(X$X[,i])[1])] <- 1
        impurete <- impurity_split(Y,split[[i]], timeScale)
        impur[i] <- impurete$impur
        toutes_imp[[i]] <- impurete$imp_list
      }
      
      if (length(unique(X$X[,i]))==1) {
        impur[i] <- Inf
        split[[i]] <- Inf
      }
    }
  }
  if (dim(X$X)[1]<=deep){
    return(list(Pure=TRUE))
  }
  if (length(unique(impur))==1 & is.element(Inf,impur)==TRUE){
    return(list(Pure=TRUE))
  }
  true_split <- which.min(impur)
  split <- split[[true_split]]
  return(list(split=split, impurete=min(impur),impur_list = toutes_imp[[true_split]], variable=which.min(impur), Pure=Pure))
}

#' Maximal Frechet tree
#'
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Shape [list]:
#' @param Image [list]:
#' @param Y [list]:
#' @param timeScale [list]:
#' @param ... : Optional parameters
#'
#' @import kmlShape
#' @import stringr
#' @import Evomorph
#' @import geomorph
#' @import survival
#'
#' @keywords internal
#'
Tmax <- function(Curve=NULL, Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL ,Y, timeScale=0.1, ...){
  
  if (is.null(Curve)==FALSE){
    Curve <- list(type="curve",X=Curve$X,id=Curve$id,time=Curve$time)
  }
  if (is.null(Scalar)==FALSE){
    Scalar <- list(type="scalar",X=Scalar$X,id=Scalar$id)
  }
  if (is.null(Factor)==FALSE){
    Factor <- list(type="factor",X=Factor$X,id=Factor$id)
  }
  
  if (is.null(Shape)==FALSE){
    Shape <- list(type="shape",X=Shape$X,id=Shape$id)
    for (k in 1:dim(Shape$X)[length(dim(Shape$X))]){
      Shape$X[,,,k] <- gpagen(Shape$X[,,,k], print.progress = FALSE)$coords
    }
  }
  
  if (is.null(Image)==FALSE){
    Image <- list(type="image",X=Image$X,id=Image$id)
  }
  # On commence par lire les arguments :
  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs
  # On va les lires en mettant la maj sur les differents elements qui le constituent :
  
  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }
  
  
  id_feuille <- rep(1,length(Y$id))
  id_feuille_prime <- id_feuille
  V_split <- NULL
  hist_nodes <- list()
  decoupe <- 1
  imp_nodes <- list()
  imp_nodes[[1]] = Inf
  impurete = Inf
  if (Y$type!="surv"){
    impurete <- impurity(Y, timeScale)
    imp_nodes[[1]] <- impurete
  }
  
  Y_pred <- list()
  Y_pred_imputation <- list()
  hist_imp_nodes <- as.matrix(cbind(1, impurete,length(unique(Y$id))))
  
  for (p in 1:(length(unique(Y$id))/2-1)){
    count_split <- 0
    for (i in 1:length(unique(id_feuille))){
      
      w_Y <- which(id_feuille==unique(id_feuille)[i])
      #On rÃ©cupÃ¨re les identifiants des individus dans ces feuilles
      which_feuilles <- unique(Y$id[w_Y])
      
      ### Il faut trouver les moyen d'avoir de manière automatique les indexes de chaque entree
      
      w_XCurve <- NULL
      w_XScalar <- NULL
      w_XFactor <- NULL
      for (l in which_feuilles){
        
        if (is.element("curve",inputs)==TRUE)  w_XCurve <- c(w_XCurve, which(Curve$id==l))
        if (is.element("scalar",inputs)==TRUE) w_XScalar <- c(w_XScalar, which(Scalar$id==l))
        if (is.element("factor",inputs)==TRUE) w_XFactor <- c(w_XFactor, which(Factor$id==l))
      }
      
      if (length(unique(Y$id[w_Y]))>1 & imp_nodes[[unique(id_feuille)[i]]] >0){
        
        #Il nous faut les entrees :
        
        if (is.element("curve",inputs)==TRUE) Curve_courant <- list(type=Curve$type, X=Curve$X[w_XCurve,],id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
        if (is.element("scalar",inputs)==TRUE) {Scalar_courant <- list(type=Scalar$type, X=Scalar$X[w_XScalar,], id=Scalar$id[w_XScalar])}
        if (is.element("factor",inputs)==TRUE) {Factor_courant <- list(type=Factor$type, X=Factor$X[w_XFactor,], id=Factor$id[w_XFactor])}
        
        #Il nous faut les sorties maintenant
        
        if (Y$type=="curve"){Y_courant <- list(type=Y$type, Y=Y$Y[w_Y], id=Y$id[w_Y], time=Y$time[w_Y])}
        if (Y$type=="factor" || Y$type=="scalar"){Y_courant <- list(type=Y$type, Y=Y$Y[w_Y], id=Y$id[w_Y])}
        if (Y$type=="shape") { Y_courant <- list(type="shape",Y=Y$Y[,,w_Y], id=Y$id[w_Y])}
        if (Y$type=="image"){Y_courant <- list(type="image",Y=Y$Y[w_Y,], id=Y$id[w_Y])}
        if (Y$type=="surv"){Y_courant <- list(type="surv", Y=Y$Y[w_Y], id=Y$id[w_Y], time=Y$time[w_Y])}
        if (Y$type=="distribution"){Y_courant <- list(type="distribution", Y=Y$Y[w_Y,], id=Y$id[w_Y])}
        if (Y$type=="sphere"){Y_courant <- list(type="sphere", Y=Y$Y[w_Y,], id=Y$id[w_Y])}
        if (Y$type=="sphere2"){Y_courant <- list(type="sphere2", Y=Y$Y[w_Y,], id=Y$id[w_Y])}
        if (Y$type=="covariance"){Y_courant <- list(type="covariance", Y=Y$Y[,,w_Y], id=Y$id[w_Y])}
        # Il nous faut maintenant faire le split sur toutes les differents types :
        
        F_SPLIT <- NULL
        decoupe <- 0
        
        if (is.element("curve",inputs)==TRUE){
          feuille_split_Curve <- var_split(Curve_courant,Y_courant, timeScale)
          F_SPLIT <- rbind(F_SPLIT,c("Curve",feuille_split_Curve$impurete))
          decoupe <- decoupe+1
        }
        
        if (is.element("scalar",inputs)==TRUE){
          feuille_split_Scalar <- var_split(Scalar_courant,Y_courant,timeScale)
          if (feuille_split_Scalar$Pure==FALSE){
            F_SPLIT <- rbind(F_SPLIT,c("Scalar",feuille_split_Scalar$impurete))
            decoupe <- decoupe +1
          }
        }
        
        if (is.element("factor",inputs)==TRUE){
          feuille_split_Factor <- var_split(Factor_courant,Y_courant,timeScale)
          if (feuille_split_Factor$Pure==FALSE){
            F_SPLIT <- rbind(F_SPLIT,c("Factor",feuille_split_Factor$impurete))
            decoupe <- decoupe +1
          }
        }
        
        
        if (decoupe >0){
          TYPE <- F_SPLIT[which.min(F_SPLIT[,2]),1]
          X <- get(TYPE)
          w_X <- get(paste("w_X",TYPE, sep=""))
          # Ensuite on va tout comparer
          
          feuille_split <- get(paste("feuille_split_",TYPE, sep="")) ####  renvoie l'impurete ainsi que le split ? gauche et a droite
          
          #if (imp_apres_split<imp_avant_split){
          
          if (Y$type=="curve"){
            Y_pred[[unique(id_feuille)[i]]] <- kmlShape::meanFrechet(as.data.frame(cbind(Y$id[w_Y], Y$time[w_Y], Y$Y[w_Y])))}
          
          if (Y$type=="scalar"){
            Y_pred[[unique(id_feuille)[i]]] <- mean(Y$Y[w_Y])
          }
          
          if (Y$type=="shape"){
            Y_pred[[unique(id_feuille)[i]]] <- mshape(Y$Y[,,w_Y, drop=FALSE])
          }
          
          if (Y$type=="factor"){ #On va renvoyer la classe dominante ::
            Table <- which.max(table(Y$Y[w_Y]))
            Y_pred[[unique(id_feuille)[i]]] <- as.factor(attributes(Table)$names)
          }
          
          if (Y$type=="image"){
            Y_pred[[unique(id_feuille)[i]]] <- apply(Y$Y[w_Y,, drop=FALSE], 2, "mean")
          }
          
          if (Y$type=="distribution"){
            Y_pred[[unique(id_feuille)[i]]] <- colMeans(Y$Y[w_Y,,drop=FALSE])
          }

          if (Y$type=="sphere"){
            S <- list()
            nRow <- nrow(Y$Y[w_Y,,drop=FALSE])
            for (h in 1:nRow){
              S[[h]] = Y$Y[w_Y,,drop=FALSE][h,] # project onto Sphere
            }
            S = riemfactory(S, name="sphere") 
            Y_pred[[unique(id_feuille)[i]]] <- as.vector(rbase.mean(S)$x)
          }

          if (Y$type=="sphere2"){
            Y_pred[[unique(id_feuille)[i]]] <- IntrinsicMean(Y$Y[w_Y,,drop=FALSE])
          }

          if (Y$type=="covariance"){
            Y_pred[[unique(id_feuille)[i]]] <- estcov(Y$Y[,,w_Y, drop=FALSE], method= method, alpha=alpha)$mean
          }
          
          if (Y$type=="surv"){
            donnees = survfit(Surv(Y$time[w_Y], Y$Y[w_Y])~1)
            Y_pred[[unique(id_feuille)[i]]] <- data.frame(times=donnees$time, traj=donnees$surv)
          }
          
          if (Y$type=="surv"){
            imp_nodes[[2*unique(id_feuille)[i]]] <- Inf
            imp_nodes[[2*unique(id_feuille)[i]+1]] <- Inf
          }
          else {
            imp_nodes[[2*unique(id_feuille)[i]]] <- feuille_split$impur_list[[1]]
            imp_nodes[[2*unique(id_feuille)[i]+1]] <- feuille_split$impur_list[[2]]
          }
          
          
          hist_imp_nodes <- rbind(hist_imp_nodes, c(2*unique(id_feuille)[i],imp_nodes[[2*unique(id_feuille)[i]]], length(which(feuille_split$split==1))))
          hist_imp_nodes <- rbind(hist_imp_nodes, c(2*unique(id_feuille)[i]+1,imp_nodes[[2*unique(id_feuille)[i]+1]], length(which(feuille_split$split==2))))
          
          
          gauche_id <- unique(Y$id[w_Y])[which(feuille_split$split==1)]
          droit_id <- unique(Y$id[w_Y])[which(feuille_split$split==2)]
          
          V_split <- rbind(V_split,c(TYPE,unique(id_feuille)[i],feuille_split$variable))
          wY_gauche <- NULL
          wY_droit <- NULL
          w_gauche <- NULL
          w_droit <- NULL
          
          for (k in 1:length(gauche_id)){
            w_gauche <- c(w_gauche, which(X$id==gauche_id[k]))
            wY_gauche <- c(wY_gauche, which(Y$id==gauche_id[k]))
          }
          
          for (k in 1:length(droit_id)){
            w_droit <- c(w_droit, which(X$id==droit_id[k]))
            wY_droit <- c(wY_droit, which(Y$id==droit_id[k]))
          }
          
          id_feuille_prime[wY_gauche] <- 2*(unique(id_feuille)[i])
          id_feuille_prime[wY_droit] <- 2*(unique(id_feuille)[i])+1
          
          if (X$type=="curve"){
            trajG <- as.data.frame(cbind(X$id[w_gauche], X$time[w_gauche], X$X[w_gauche,feuille_split$variable]))
            trajD <- as.data.frame(cbind(X$id[w_droit], X$time[w_droit], X$X[w_droit,feuille_split$variable]))
            meanFg <- as.matrix(kmlShape::meanFrechet(trajG))
            meanFd <- as.matrix(kmlShape::meanFrechet(trajD))
          }
          
          if (X$type=="factor"){
            meanFg <- unique(X$X[w_gauche, feuille_split$variable])
            meanFd <- unique(X$X[w_droit,feuille_split$variable])
          }
          
          if (X$type=="scalar"){
            meanFg <- mean(X$X[w_gauche,feuille_split$variable])
            meanFd <- mean(X$X[w_droit,feuille_split$variable])
          }
          
          hist_nodes[[2*(unique(id_feuille)[i])]] <- meanFg
          hist_nodes[[2*(unique(id_feuille)[i])+1]] <- meanFd
          
          count_split <- count_split+1
          
          # QUelles sont les feuilles courantes
          feuilles_courantes <- unique(id_feuille_prime)
          info_feuilles <- hist_imp_nodes[is.element(hist_imp_nodes[,1], feuilles_courantes),]
          impurete <- c(impurete, sum(info_feuilles[,2]*info_feuilles[,3]/length(unique(Y$id))))
        }
        
      }
    }
    
    id_feuille <- id_feuille_prime
    
    if (count_split ==0 ){
      V_split <- data.frame(V_split)
      names(V_split) <- c("type","num_noeud", "var_split")
      for (q in unique(id_feuille)){
        w <- which(id_feuille == q)
        if (Y$type=="curve"){
          Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y$id[w], Y$time[w], Y$Y[w]))
        }
        if(Y$type=="scalar"){
          Y_pred[[q]]<- mean(Y$Y[w])
        }
        if(Y$type=="factor"){
          Table <- which.max(table(Y$Y[w]))
          Y_pred[[q]] <-  as.factor(attributes(Table)$names)
        }
        
        if (Y$type=="shape"){
          Y_pred[[q]] <- mshape(Y$Y[,,w, drop=FALSE])
        }
        
        if (Y$type=="image"){
          Y_pred[[q]] <- apply(Y$Y[w,, drop=FALSE], 2, "mean")
        }
        
        if (Y$type=="distribution"){
          Y_pred[[q]] <- colMeans(Y$Y[w,,drop=FALSE])
        }

        if (Y$type=="sphere"){
          S <- list()
          nRow <- nrow(Y$Y[w,,drop=FALSE])
          for (h in 1:nRow){
            S[[h]] = Y$Y[w,,drop=FALSE][h,] # project onto Sphere
          }
          S = riemfactory(S, name="sphere")  
          Y_pred[[q]] <- as.vector(rbase.mean(S)$x)
        }

        if (Y$type=="sphere2"){
          Y_pred[[q]] <- IntrinsicMean(Y$Y[w,,drop=FALSE])
        }

        if (Y$type=="covariance"){
          Y_pred[[q]] <- estcov(Y$Y[,,w,drop=FALSE], method= method, alpha=alpha)$mean
        }
        
      }
      if (Y$type=="factor"){
        Ylevels <- unique(Y$Y)
        return(list(feuilles = id_feuille,  V_split=V_split, impurity=impurete, hist_nodes=hist_nodes, Y_pred = Y_pred, Y=Y, hist_imp_nodes=hist_imp_nodes, Alpha =0, Ylevels=Ylevels))
      }
      return(list(feuilles = id_feuille,  V_split=V_split, impurity=impurete, hist_nodes=hist_nodes, Y_pred = Y_pred, Y=Y, hist_imp_nodes=hist_imp_nodes, Alpha =0))
    }
  }
  
  V_split <- data.frame(V_split)
  names(V_split) <- c("type","num_noeud", "var_split")
  for (q in unique(id_feuille)){
    w <- which(id_feuille == q)
    if (Y$type=="curve"){
      Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y$id[w], Y$time[w], Y$Y[w]))
    }
    if(Y$type=="scalar"){
      Y_pred[[q]]<- mean(Y$Y[w])
    }
    if(Y$type=="factor"){
      Table <- which.max(table(Y$Y[w]))
      Y_pred[[q]] <-  as.factor(attributes(Table)$names)
    }
    if (Y$type=="shape"){
      Y_pred[[q]] <- mshape(Y$Y[,,w, drop=FALSE])
    }
    
    if (Y$type=="image"){
      Y_pred[[q]] <- apply(Y$Y[w,, drop=FALSE], 2, "mean")
    }
    
    if (Y$type=="surv"){
      donnees = survfit(Surv(Y$time[w], Y$Y[w])~1)
      Y_pred[[q]] <- data.frame(times=donnees$time, traj=donnees$surv)
    }
    
    if (Y$type=="distribution"){
      Y_pred[[q]] <- colMeans(Y$Y[w,,drop=FALSE])
    }

    if (Y$type=="sphere"){
      S <- list()
      nRow <- nrow(Y$Y[w,,drop=FALSE])
      for (h in 1:nRow){
        S[[h]] = Y$Y[w,,drop=FALSE][h,] # project onto Sphere
      }
      S = riemfactory(S, name="sphere")  
      Y_pred[[q]] <- as.vector(rbase.mean(S)$x)
    }

    if (Y$type=="sphere2"){
      Y_pred[[q]] <- IntrinsicMean(Y$Y[w,,drop=FALSE])
    }

    if (Y$type=="covariance"){ 
      Y_pred[[q]] <- estcov(Y$Y[,,w, drop=FALSE], method= method, alpha=alpha)$mean
    }
    
  }
  if (Y$type=="factor"){
    Ylevels <- unique(Y$Y)
    return(list(feuilles = id_feuille, V_split=V_split, impurity=impurete, hist_nodes=hist_nodes, Y_pred = Y_pred, Y=Y, hist_imp_nodes=hist_imp_nodes, Alpha =0, Ylevels=Ylevels))
  }
  return(list(feuilles = id_feuille, V_split=V_split, impuity=impurete, hist_nodes=hist_nodes, Y_pred= Y_pred, Y=Y, hist_imp_nodes=hist_imp_nodes, Alpha=0))
}

#' Sub trees  extractor
#'
#' @param tree
#' @param t
#'
#'
#' @keywords internal
branche <- function(tree, t){
  Y <- list()
  f <- unique(tree$feuilles)
  sous_split <- tree$V_split[which(tree$V_split[,2]==t),]
  N <- 2
  g <- which(tree$V_split[,2]==2*t)
  d <- which(tree$V_split[,2]==2*t+1)
  noeuds_courants <- as.numeric(as.character(tree$V_split[c(g,d),2]))
  noeuds_courants1 <- noeuds_courants
  sous_split <- rbind(sous_split, tree$V_split[c(g,d),])
  sous_feuilles <- NULL
  hist_nodes <- list()
  if (length(g)>0) {hist_nodes[[2*t]] <- tree$hist_nodes[[2*t]]}
  if (length(d)>0) {hist_nodes[[2*t+1]] <- tree$hist_nodes[[2*t+1]]}
  if (length(d)== 0) {sous_feuilles <- c(sous_feuilles, 2*t+1)
  Y[[2*t+1]] <- tree$Y_pred[[2*t+1]]}
  if (length(g)== 0) {sous_feuilles <- c(sous_feuilles, 2*t)
  Y[[2*t]] <- tree$Y_pred[[2*t]]}
  racine <- t
  if (length(noeuds_courants)>0) {
    while(N>0){
      p <- 0
      courant_prime <- NULL
      for (l in noeuds_courants){
        g <- which(tree$V_split[,2]==2*l)
        d <- which(tree$V_split[,2]==2*l+1)
        
        if (length(g)>0){ p <- p+2
        courant_prime <- c(courant_prime, as.numeric(as.character(tree$V_split[g,2])))
        sous_split <- rbind(sous_split, tree$V_split[g,])
        hist_nodes[[2*l]] <- tree$hist_nodes[[2*l]]}
        
        if (length(d)>0){ p <- p+2
        courant_prime <- c(courant_prime, as.numeric(as.character(tree$V_split[d,2])))
        sous_split <- rbind(sous_split, tree$V_split[d,])
        hist_nodes[[2*l+1]] <- tree$hist_nodes[[2*l+1]]}
        
        if(length(g)==0) {sous_feuilles <- c(sous_feuilles,2*l)
        Y[[2*l]] <- tree$Y_pred[[2*l]]}
        
        if (length(d)==0) { sous_feuilles <- c(sous_feuilles, 2*l+1)
        Y[[2*l+1]] <- tree$Y_pred[[2*l+1]]}
      }
      noeuds_courants <- courant_prime
      N <-p
    }
  }
  
  if (length(noeuds_courants1)==0) {sous_feuilles <- c(2*t, 2*t+1)}
  
  ## C'est maintenant que ca devient coton :::
  # Il faut recuperer les id des gens qui sont
  
  s_feuilles <- NULL
  s_id <- NULL
  s_time <- NULL
  s_Y <- NULL
  
  for(f in unique(sous_feuilles)){
    w <- which(tree$feuilles==f)
    s_feuilles <- c(s_feuilles, tree$feuilles[w])
    s_id <- c(s_id, tree$Y$id[w])
    if (tree$Y$type=="curve"){
      s_time <- c(s_time,tree$Y$time[w])
    }
    #s_time <- c(s_time, tree$time[w])
    if (tree$Y$type=="shape" || tree$Y$type=="image" || tree$Y$type=="distribution" || tree$Y$type=="covariance" || tree$Y$type=="sphere" || tree$Y$type=="sphere2") s_Y <- c(s_Y,w)
    else s_Y <- c(s_Y, tree$Y$Y[w])
  }
  if (tree$Y$type=="shape" || tree$Y$type=="image" || tree$Y$type=="covariance") s_Y <- tree$Y$Y[,,s_Y,drop=FALSE]
  if (tree$Y$type=="distribution" || tree$Y$type=="sphere" || tree$Y$type=="sphere2") s_Y <- tree$Y$Y[s_Y,,drop=FALSE]
  #### il faut maintenant calculer l'impurete de la branche ainsi que celle du noeud t
  #### impurete dans le noeud racine :::
  impurity_racine <- tree$hist_imp_nodes[which(tree$hist_imp_nodes[,1]==racine),2]
  n_racine <- tree$hist_imp_nodes[which(tree$hist_imp_nodes[,1]==racine),3]
  n_base <- tree$hist_imp_nodes[1,3]
  impurity_racine <- impurity_racine*(n_racine/n_base)
  
  impurity_T <- 0
  for (i in unique(s_feuilles)){
    w <- which(tree$hist_imp_nodes[,1]==i)
    prop <- tree$hist_imp_nodes[w,3]/n_base
    impurity_T <- impurity_T + tree$hist_imp_nodes[w,2]*prop
  }
  if (tree$Y$type=="curve"){
    sous_Y <- list(type=tree$Y$type, Y=s_Y, id = s_id, time=s_time)
  }
  else sous_Y <- list(type=tree$Y$type, Y=s_Y, id = s_id)
  return(list(feuilles=s_feuilles, V_split = sous_split, hist_nodes=hist_nodes, Y=sous_Y, impurity_T = impurity_T, impurity_racine = impurity_racine, n_racine=n_racine, Y_pred=Y))
}

#' Detect and destroy nodes
#'
#' @param tree
#'
#'
#' @keywords internal
noeuds_deg <- function(tree){
  noeuds <- as.numeric(as.character(tree$V_split$num_noeud))
  deg <- NULL
  alpha <- rep()
  mat_pen <- matrix(0, length(noeuds), 5)
  mat_pen[,1] <- noeuds
  for (t in noeuds){
    b <- branche(tree,t) ### on recupère la branche associee à t
    if (length(unique(b$feuilles))>1){
      mat_pen[which(noeuds==t), 2] <- b$impurity_racine
      mat_pen[which(noeuds==t), 3] <- b$impurity_T
      mat_pen[which(noeuds==t), 4] <- length(unique(b$feuilles))
      mat_pen[which(noeuds==t), 5] <- (b$impurity_racine-b$impurity_T)/(length(unique(b$feuilles))-1)}
    #pen <- mat_pen[which(noeuds==t), 5]
    #err <- b$impurity_T + pen*length(unique(b$feuilles)) - b$impurity_racine - pen
    #print(err)
  }
  alpha <- min(mat_pen[,5])
  err <- rep(0, length(noeuds))
  for (i in  1:dim(mat_pen)[1]){
    err[i] <- round(mat_pen[i,3] + alpha*mat_pen[i,4] - mat_pen[i,2] - alpha, 5)
    if (err[i]==0){
      deg <- rbind(deg, c(mat_pen[i,1], alpha))
    }
  }
  return(deg)
}

#' General pruning function
#'
#' @param tree
#'
#'
#' @keywords internal
elagage <- function(tree){
  
  t_feuilles <- NULL
  t_id <- NULL
  t_time <- NULL
  t_split <- NULL
  t_hist <- NULL
  t_Y <- tree$Y
  
  
  tree_courant <- tree
  nb_feuilles <- length(unique(tree$feuilles))
  n_max <- nb_feuilles
  courant <- 2
  TREES <- list()
  TREES[[1]] <- tree
  ##### il faut aussi trouver les d?coupe superficielles :::: on garde un historique des d?coupes :::::
  while(nb_feuilles >1){
    deg <- noeuds_deg(tree_courant)
    if (dim(deg)[1]>1) deg <- apply(deg, 2, sort, decreasing=TRUE)
    t_feuilles_courant <- tree_courant$feuilles
    for (t in deg[,1]){
      b <- branche(tree_courant, t)
      feuilles_b <- unique(b$feuilles)
      w_feuilles <- NULL
      for (f in feuilles_b ){
        w_feuilles <- c(w_feuilles, which(tree_courant$feuilles==f))
      }
      t_feuilles_courant[w_feuilles] <- t
      #### il faut maintenant retirer toute la branche de t
      
      nodes <- as.numeric(as.character(b$V_split[,2]))
      w_nodes <- NULL
      for (node in nodes){
        w_nodes <- c(w_nodes, which(tree_courant$V_split[,2]==node))
      }
      
      t_split_courant <- tree_courant$V_split[-w_nodes,, drop = FALSE]
      ##### il faut alors recalculer l'importance dans les feuilles
      tree_courant <- list(feuilles=t_feuilles_courant, V_split = t_split_courant,hist_nodes=tree$hist_nodes, Y=tree$Y, hist_imp_nodes=tree$hist_imp_nodes, Y_pred = tree$Y_pred, Alpha=unique(deg[,2]))
      TREES[[courant]] <- tree_courant
    }
    courant <- courant+1
    nb_feuilles <- length(unique(tree_courant$feuilles))
  }
  
  return(TREES)
}

#' General Frechet Tree
#'
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Y [list]:
#' @param timeScale [numeric]:
#' @param ncores [numeric]:
#'
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import kmlShape
#' @import Evomorph
#' @import pbapply
#'
#' @export
#'
FrechetTree <- function(Curve=NULL,Scalar=NULL,Factor=NULL,Y,timeScale=0.1, ncores=NULL){
  
  ### Il faut normaliser les elements des formes ::
  
  if (is.null(ncores)==TRUE) ncores <- detectCores()-1
  
  if (Y$type=="shape"){
    Y$Y <- gpagen(Y$Y,print.progress = FALSE)$coords
  }
  
  TMAX <- Tmax(Curve=Curve,Scalar = Scalar,Factor=Factor, Shape=NULL,Image=NULL,Y,timeScale = timeScale)
  
  if (Y$type=="shape") dime <- dim(Y$Y)[1:2]
  
  
  elag_max <- elagage(TMAX)
  ALPHA <- rep(NA, length(elag_max))
  for (i in 1:length(ALPHA)){
    ALPHA[i] <- elag_max[[i]]$Alpha
  }
  #### on transforme le tout en beta
  beta <- rep(NA, length(ALPHA))
  beta[length(ALPHA)] <- ALPHA[length(ALPHA)]
  for (i in 1:(length(ALPHA)-1)){
    beta[i] <- sqrt(abs(ALPHA[i]*ALPHA[i+1]))
  }
  
  #### Il faut faire faire les sous ensemble de validation crois?e::::
  #ELAG <- list()
  n_folds <- 10
  VC <- sample(rep(1:n_folds, length.out= length(unique(Y$id))))
  #tmax <- list()
  #APP <- list()
  #err <- matrix(0, length(beta), n_folds)
  
  
  Scalar.app <- Scalar
  Curve.app <- Curve
  Factor.app <- Factor
  
  Scalar.val <- NULL
  Factor.val <- NULL
  Curve.val <- NULL
  
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  #p=1
  #err <- foreach(p = 1:n_folds, .packages='kmlShape', .combine="cbind") %dopar% {
  
  err <- pbsapply(1:n_folds, FUN=function(i){
    
    
    res <- rep(NA,length(beta))
    app <- unique(Y$id)[which(VC!=i)] ### on r?cup?re les identifiants
    w <- NULL
    wCurve <- NULL
    wScalar <- NULL
    wFactor <- NULL
    
    for (a in app){
      w <- c(w, which(Y$id==a))
      if (is.null(Scalar)!=TRUE) wScalar <- c(wScalar, which(Scalar$id==a))
      if (is.null(Factor)!=TRUE) wFactor <- c(wFactor, which(Factor$id==a))
      if (is.null(Curve)!=TRUE) wCurve <- c(wCurve, which(Curve$id==a))
      
    }
    APP <- w
    
    ### On prend les elements d'apprentissage maintenant :::
    
    if (is.null(Scalar)!=TRUE){
      Scalar.app <- list(X=Scalar$X[wScalar,,drop=FALSE], id=Scalar$id[wScalar])
      Scalar.val <- list(type="scalar",X=Scalar$X[-wScalar,,drop=FALSE], id=Scalar$id[-wScalar])
    }
    
    
    if (is.null(Factor)!=TRUE){
      Factor.app <- list(X=Factor$X[wFactor,,drop=FALSE], id=Factor$id[wFactor])
      Factor.val <- list(type="factor",X=Factor$X[-wFactor,,drop=FALSE], id=Factor$id[-wFactor])
    }
    
    if (is.null(Curve)!=TRUE){
      Curve.app <- list(X=Curve$X[wCurve,,drop=FALSE], id=Curve$id[wCurve], time=Curve$time[wCurve])
      Curve.val <- list(type="curve",X=Curve$X[-wCurve,,drop=FALSE], id=Curve$id[-wCurve], time=Curve$time[-wCurve])
    }
    
    if (Y$type=="curve"){
      Y.app <- list(type="curve",Y=Y$Y[w],id=Y$id[w],time=Y$time[w])
      Y.val <- list(type="curve",Y=Y$Y[-w],id=Y$id[-w],time=Y$time[-w])
    }
    
    if (Y$type=="factor" || Y$type=="scalar"){
      Y.app <- list(type=Y$type,Y=Y$Y[w],id=Y$id[w],time=Y$time[w])
      Y.val <- list(type=Y$type,Y=Y$Y[-w],id=Y$id[-w],time=Y$time[-w])
    }
    
    if (Y$type=="shape"){
      Y.app <- list(type=Y$type,Y=Y$Y[,,w],id=Y$id[w])
      Y.val <- list(type=Y$type,Y=Y$Y[,, -w],id=Y$id[-w])
    }
    
    if (Y$type=="image"){
      Y.app <- list(type=Y$type,Y=Y$Y[w,],id=Y$id[w])
      Y.val <- list(type=Y$type,Y=Y$Y[-w,],id=Y$id[-w])
    }
    
    if (Y$type=="distribution"){
      Y.app <- list(type=Y$type,Y=Y$Y[w,],id=Y$id[w])
      Y.val <- list(type=Y$type,Y=Y$Y[-w,],id=Y$id[-w])
      if(is.vector(Y.val$Y)) Y.val$Y=matrix(Y.val$Y,nrow = 1)
    }

    if (Y$type=="sphere"){
      Y.app <- list(type=Y$type,Y=Y$Y[w,],id=Y$id[w])
      Y.val <- list(type=Y$type,Y=Y$Y[-w,],id=Y$id[-w])
      if(is.vector(Y.val$Y)) Y.val$Y=matrix(Y.val$Y,nrow = 1)
    }

    if (Y$type=="sphere2"){
      Y.app <- list(type=Y$type,Y=Y$Y[w,],id=Y$id[w])
      Y.val <- list(type=Y$type,Y=Y$Y[-w,],id=Y$id[-w])
      if(is.vector(Y.val$Y)) Y.val$Y=matrix(Y.val$Y,nrow = 1)
    }

    if (Y$type=="covariance"){
      Y.app <- list(type=Y$type,Y=Y$Y[,,w],id=Y$id[w])
      Y.val <- list(type=Y$type,Y=Y$Y[,,-w],id=Y$id[-w])
      if(is.matrix(Y.val$Y)) Y.val$Y=array(Y.val$Y, dim = c(dim(Y.val$Y)[1],dim(Y.val$Y)[2],1))
    }
    
    tmax <- Tmax(Curve = Curve.app,Scalar = Scalar.app,Factor = Factor.app,Y=Y.app, timeScale = timeScale)
    
    ELAG <- elagage(tmax)
    
    pen <- rep(NA,length(ELAG))
    #pen[length(ELAG[[p]])] <- ELAG[[p]][[length(ELAG[[p]])]]$Alpha
    for (l in 1:length(pen)){
      pen[l] <- ELAG[[l]]$Alpha
    }
    
    for (k in 1:length(beta)){
      sous_arbre <- ELAG[[which.min(abs(pen-beta[k]))]]
      where <- pred.FT(sous_arbre,Curve = Curve.val,Scalar=Scalar.val,Factor = Factor.val,timeScale = timeScale) #### on doit trouver les feuilles de pr?diction :::
      ##### il nous faut maintenant pr?dire les diff?rentes courbes ::::
      err_courant <- rep(0, length(where))
      
      for (j in 1:length(where)){
        ww <- which(Y.val$id == unique(Y.val$id)[j])
        #mean_courant <- DouglasPeuckerEpsilon(sous_arbre$Y_curves[[where[j]]][,1],sous_arbre$Y_curves[[where[j]]][,2], 0.01)
        ## On regarde si on a bien une sortie qui est une courbe:
        if (Y$type=="curve") err_courant[j] <-  kmlShape::distFrechet(Y.val$time[ww], Y.val$Y[ww],sous_arbre$Y_pred[[where[j]]][,1], sous_arbre$Y_pred[[where[j]]][,2])
        if (Y$type=="scalar") err_courant[j] <- (Y.val$Y[ww]-where[j])^2
        if (Y$type=="factor") err_courant[j] <- 1*(Y.val$Y[ww]==where[j])
        if (Y$type=="image") err_courant[j] <- mean((Y.val[ww,]-sous_arbre$Y_pred[[where[j]]])^2)
        if (Y$type=="shape") err_courant[j] <- ShapeDist(array(sous_arbre$Y_pred[[where[j]]],dim=c(dime[1],dime[2],1)),Y.val$Y[,,ww])
        if (Y$type=="distribution") {
          #print(Y.val$Y); print(sous_arbre$Y_pred[[where[j]]])
          err_courant[j] <- sqrt(pracma::trapz(qSup, (Y.val$Y[ww,] - sous_arbre$Y_pred[[where[j]]])^2))
        }
        if (Y$type=="sphere") {
          #print(Y.val$Y); print(sous_arbre$Y_pred[[where[j]]])
          tmp = crossprod(Y.val$Y[ww,],sous_arbre$Y_pred[[where[j]]])[1,1]
          tmp = sign(tmp) * min(1, abs(tmp))
          err_courant[j] <- acos(tmp)
        }
        if (Y$type=="sphere2") {
          #print(Y.val$Y); print(sous_arbre$Y_pred[[where[j]]])
          err_courant[j] <- distCosine(Y.val$Y[ww,],sous_arbre$Y_pred[[where[j]]],r=1)
        }
        if (Y$type=="covariance") {
          #print(Y.val$Y); print(sous_arbre$Y_pred[[where[j]]])
          err_courant[j] <- distcov(Y.val$Y[,,ww], sous_arbre$Y_pred[[where[j]]], method=method, alpha=alpha)
        }
      }                                               
      res[k] <- mean(err_courant)
    }
    res
  },cl=cl)
  
  parallel::stopCluster(cl)
  
  SD <- apply(err, 1, "sd")
  err_M <- apply(err, 1, "mean")
  
  ### On prend le meilleur modèle:
  seuil <- min(err_M) + SD[which.min(err_M)]
  ### il faut s?lectionner le meilleur arbre
  optimal.tree <- max(which(err_M<=seuil))
  
  beta.opt <- beta[optimal.tree]
  final_tree <- elag_max[[optimal.tree]]
  
  #### On va s?lectionner l'arbre optimal pour chaque ensemble d'apprentissage puis calculer l'importance des variables sur ceux-ci::
  #Importance <- matrix(0, n_folds, dim(X)[2])
  #err_arbres_select <- rep(NA, n_folds)
  #for (k in 1:n_folds){
  ### on r?cup?re les ?l?ments de validation::::
  #X.val <- X[-APP[[k]],]
  #Y.val <- Y[-APP[[k]]]
  #time.val <- time[-APP[[k]]]
  #id.val <- id[-APP[[k]]]
  
  # pen <- rep(NA,length(ELAG[[k]]))
  #for (l in 1:length(pen)){
  #  pen[l] <- ELAG[[k]][[l]]$Alpha
  #}
  #}
  ## On va faire l'affichage de la selection de l'abre
  #plot(err_M)
  #lines(rep(seuil, length(err_M)),col=2)
  #points(optimal.tree,err_M[optimal.tree], col=2)
  
  m_leafs <- max(unique(final_tree$feuilles))
  return(list(feuilles = final_tree$feuilles, V_split=final_tree$V_split, hist_nodes=final_tree$hist_nodes[1:m_leafs], Y_pred=final_tree$Y_pred[1:m_leafs], err_elag = err_M, seuil=seuil, Y=Y))
}

#' Predict Frechet tree
#'
#' @param tree : Frechet tree.
#' @param Curve [list]: A list that contains the input curves.
#' @param Scalar [list]: A list that contains the input scalars.
#' @param Factor [list]: A list that contains the input factors.
#' @param Shape [list]: A list that contains the input shape.
#' @param Image [list]: A list that contains the input images.
#' @param aligned.shape [logical]: \code{TRUE} if the input shapes are aligned and normalized (\code{aligned.shape=FALSE} by default)
#' @param timeScale [numeric]: Time scale for the input and output curves (\code{timeScale=0.1} by default)
#'
#' @import stringr
#' @import geomorph
#' @import kmlShape
#' @import Evomorph
#'
#' @export
#'
pred.FT <- function(tree, Curve=NULL,Scalar=NULL,Factor=NULL,Shape=NULL,Image=NULL, aligned.shape=FALSE ,timeScale=0.1){
  
  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs
  
  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }
  
  id.pred <- unique(get(Inputs[1])$id)
  
  if (is.element("shape",inputs)==TRUE & aligned.shape==FALSE){
    for (k in 1:dim(Shape$X)[length(dim(Shape$X))]){
      Shape$X[,,,k] <- gpagen(Shape$X[,,,k],print.progress = FALSE)$coords
    }
  }
  
  
  if (tree$Y$type=="factor"){
    pred <- factor(rep(NA, length(id.pred)),levels=tree$Ylevels)
  }
  
  else{pred <- rep(NA,length(id.pred))}
  
  for (i in 1:length(id.pred)){
    
    if (is.element("curve",inputs)==TRUE) wCurve <- which(Curve$id==id.pred[i])
    if (is.element("scalar",inputs)==TRUE) wScalar <- which(Scalar$id==id.pred[i])
    if (is.element("factor",inputs)==TRUE) wFactor <- which(Factor$id==id.pred[i])
    if (is.element("shape",inputs)==TRUE) wShape <- which(Shape$id==id.pred[i])
    if (is.element("image",inputs)==TRUE) wImage <- which(Image$id==id.pred[i])
    
    noeud_courant <- 1
    
    while (is.element(noeud_courant, tree$feuilles)==FALSE){
      
      X <- get(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),1]))
      type <- str_to_lower(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),1]))
      var.split <- as.numeric(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),3]))
      
      # Maintenant il nous faut regarder la difference entre la moyenne à gauche et a droite et conclure :
      
      meanG <- tree$hist_nodes[[2*noeud_courant]]
      meanD <- tree$hist_nodes[[2*noeud_courant+1]]
      
      if (type=="curve"){
        distG <- distFrechet(meanG[,1], meanG[,2], X$time[wCurve], X$X[wCurve,var.split], timeScale = timeScale)
        distD <- distFrechet(meanD[,1], meanD[,2], X$time[wCurve], X$X[wCurve,var.split], timeScale = timeScale)
      }
      if (type=="scalar"){
        distG <- abs(meanG- X$X[wScalar,var.split])
        distD <- abs(meanD-X$X[wScalar,var.split])
      }
      
      if (type=="shape"){
        elementz <- array(X$X[,,wShape,var.split],dim = c(nrow(meanG),ncol(meanG),1))
        distG <- ShapeDist(elementz,meanG)
        distD <- ShapeDist(elementz, meanD)
      }
      
      if (type=="image"){
        distG <- mean((X$X[wImage,,var.split]-meanG)^2)
        distD <- mean((X$X[wImage,,var.split]-meanD)^2)
      }
      
      if (type=="factor"){
        distG <- -1*(is.element(X$X[wFactor,var.split],meanG))
        distD <- -1*(is.element(X$X[wFactor,var.split],meanD))
      }
      
      if (distG <= distD) { noeud_courant <- 2*noeud_courant}
      if (distD < distG) {noeud_courant <- 2*noeud_courant +1}
      
      
    }
    
    if(tree$Y$type=="curve" || tree$Y$type=="image" || tree$Y$type=="shape" || tree$Y$type=="surv" || tree$Y$type=="distribution" || tree$Y$type=="covariance" || tree$Y$type=="sphere" || tree$Y$type=="sphere2"){
      pred[i] <- noeud_courant
    }
    
    else{
      pred[i] <- tree$Y_pred[[noeud_courant]]
    }
  }
  return(pred)
}

#' Randomized Frechet tree
#'
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Shape [list]:
#' @param Image [list]:
#' @param Y [list]:
#' @param mtry [integer]:
#' @param ERT [logical]:
#' @param aligned.shape [logical]:
#' @param timeScale [numeric]:
#' @param ntry [numeric]:
#' @param ... : option
#'
#' @import kmlShape
#' @import stringr
#' @import Evomorph
#' @import geomorph
#' @import survival
#'
#' @keywords internal
Rtmax <- function(Curve=NULL, Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL,Y,mtry,ERT=FALSE,aligned.shape=FALSE,ntry=3, timeScale=0.1, ...){
  
  
  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs
  
  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }
  
  impurity_feuilles <- NULL
  V_split <- NULL
  hist_nodes <- list()
  id_boot <- unique(sample(unique(Y$id), length(unique(Y$id)), replace=TRUE))
  boot <- id_boot
  decoupe <- 1
  
  wXCurve <- NULL
  wXScalar <- NULL
  wXFactor <- NULL
  wXShape <- NULL
  wXImage <- NULL
  wY <- NULL
  
  for (k in id_boot){
    wY <- c(wY, which(Y$id==k))
    if (is.element("curve",inputs)==TRUE) wXCurve <- c(wXCurve, which(Curve$id==k))
    if (is.element("scalar",inputs)==TRUE) wXScalar <- c(wXScalar, which(Scalar$id==k))
    if (is.element("factor",inputs)==TRUE) wXFactor <- c(wXFactor, which(Factor$id==k))
    if (is.element("shape",inputs)==TRUE) wXShape <- c(wXShape, which(Shape$id==k))
    if (is.element("image",inputs)==TRUE) wXImage <- c(wXImage, which(Image$id==k))
  }
  
  Y_pred <- list()
  Y_pred_surv  <- list()
  
  if (is.element("curve",inputs)==TRUE) Curve_boot <- list(type=Curve$type,   X=Curve$X[wXCurve,, drop=FALSE], id= Curve$id[wXCurve], time = Curve$time[wXCurve]) ### bootstrap pour les courbes
  if (is.element("scalar",inputs)==TRUE) Scalar_boot <- list(type=Scalar$type,   X=Scalar$X[wXScalar,, drop=FALSE], id= Scalar$id[wXScalar]) ### bootstrap pour les courbes
  if (is.element("factor",inputs)==TRUE) Factor_boot <- list(type=Factor$type,   X=Factor$X[wXFactor,, drop=FALSE], id= Factor$id[wXFactor])
  if (is.element("shape",inputs)==TRUE) Shape_boot <- list(type=Shape$type,   X=Shape$X[,,wXShape, , drop=FALSE], id= Shape$id[wXShape])
  if (is.element("image",inputs)==TRUE) Image_boot <- list(type=Image$type,   X=Image$X[wXImage,,, drop=FALSE], id= Image$id[wXImage])
  
  
  if (Y$type=="curve" || Y$type=="surv") {Y_boot <- list(type=Y$type,Y=Y$Y[wY], id=Y$id[wY], time=Y$time[wY])} ### idem pour Y
  if (Y$type=="shape") {Y_boot <- list(type=Y$type, Y=Y$Y[,,wY], id=Y$id[wY])}
  if (Y$type=="image") {Y_boot <- list(type=Y$type, Y=Y$Y[wY,], id=Y$id[wY])}
  if (Y$type=="factor" || Y$type=="scalar") {Y_boot <- list(type=Y$type,Y=Y$Y[wY], id=Y$id[wY])}
  if (Y$type=="distribution") {Y_boot <- list(type=Y$type,Y=Y$Y[wY,], id=Y$id[wY])}
  if (Y$type=="sphere") {Y_boot <- list(type=Y$type,Y=Y$Y[wY,], id=Y$id[wY])}
  if (Y$type=="sphere2") {Y_boot <- list(type=Y$type,Y=Y$Y[wY,], id=Y$id[wY])}
  if (Y$type=="covariance") {Y_boot <- list(type=Y$type,Y=Y$Y[,,wY], id=Y$id[wY])}
  
  imp_nodes <- list()
  imp_nodes[[1]] = Inf
  impurete = Inf
  if (Y$type!="surv"){
    impurete <- impurity(Y_boot, timeScale)
    imp_nodes[[1]] <- impurete
  }
  
  id_feuille <- rep(1,length(Y_boot$id)) #### localisation des feuilles de l'arbre
  id_feuille_prime <- id_feuille
  
  for (p in 1:(length(unique(Y_boot$id))/2-1)){
    count_split <- 0
    for (i in 1:length(unique(id_feuille))){
      # Il faut que l'on regarde le tirage des variables de manière aleatoire :
      V <- NULL
      for (v in Inputs){
        V <- c(V, rep(get(v)$type,dim(get(v)$X)[length(dim(get(v)$X))]))
      }
      variables <- sample(V,mtry) # Maintenant on sait combien on doit en tirer dans chaque espace
      # On ne va regarder que les espaces tires :
      split.spaces <- unique(variables)
      
      # variables <- sample(c(1:dim(X_boot$X[,,drop=FALSE])[2]),mtry)
      w <- which(id_feuille==unique(id_feuille)[i])
      wXCurve <- NULL
      wXScalar <- NULL
      wXFactor <- NULL
      wXShape <- NULL
      wXImage <- NULL
      
      for (l in unique(Y_boot$id[w])){
        if (is.element("curve",inputs)==TRUE) wXCurve <- c(wXCurve, which(Curve_boot$id==l))
        if (is.element("scalar",inputs)==TRUE) wXScalar <- c(wXScalar, which(Scalar_boot$id==l))
        if (is.element("factor",inputs)==TRUE) wXFactor <- c(wXFactor, which(Factor_boot$id==l))
        if (is.element("shape",inputs)==TRUE) wXShape <- c(wXShape, which(Shape_boot$id==l))
        if (is.element("image",inputs)==TRUE) wXImage <- c(wXImage, which(Image_boot$id==l))
      }
      
      if (length(unique(Y_boot$id[w]))>1 & imp_nodes[[unique(id_feuille)[i]]] >0){
        
        # On est ici
        
        if (is.element("curve",split.spaces)==TRUE){
          tirageCurve <- sample(1:ncol(Curve$X),length(which(variables=="curve")))
          Curve_courant <- list(type = Curve_boot$type, X=Curve_boot$X[wXCurve,tirageCurve, drop=FALSE], id=Curve_boot$id[wXCurve, drop=FALSE], time=Curve_boot$time[wXCurve, drop=FALSE])
        }
        
        if (is.element("scalar",split.spaces)==TRUE){
          tirageScalar <- sample(1:ncol(Scalar$X),length(which(variables=="scalar")))
          Scalar_courant <- list(type = Scalar_boot$type, X=Scalar_boot$X[wXScalar,tirageScalar, drop=FALSE], id=Scalar_boot$id[wXScalar, drop=FALSE])
        }
        
        if (is.element("factor",split.spaces)==TRUE){
          tirageFactor <- sample(1:ncol(Factor$X),length(which(variables=="factor")))
          Factor_courant <- list(type = Factor_boot$type, X=Factor_boot$X[wXFactor,tirageFactor, drop=FALSE], id=Factor_boot$id[wXFactor, drop=FALSE])
        }
        
        if (is.element("shape",split.spaces)==TRUE){
          tirageShape <- sample(1:dim(Shape$X)[length(dim(Shape$X))],length(which(variables=="shape")))
          Shape_courant <- list(type = Shape_boot$type, X=Shape_boot$X[,,wXShape,tirageShape, drop=FALSE], id=Shape_boot$id[wXShape, drop=FALSE])
        }
        
        if (is.element("image",split.spaces)==TRUE){
          tirageImage <- sample(1:dim(Image$X)[3],length(which(variables=="image")))
          Image_courant <- list(type = Image_boot$type, X=Image_boot$X[wXImage,,tirageImage, drop=FALSE], id=Image_boot$id[wXImage])
        }
        
        if (Y_boot$type=="curve" || Y_boot$type=="surv"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w], id=Y_boot$id[w], time=Y_boot$time[w])
        }
        
        if (Y_boot$type=="shape"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[,,w, drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }
        
        if (Y_boot$type=="image"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w, ,drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }
        
        
        if (Y_boot$type=="factor" || Y_boot$type=="scalar"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w, drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }
        
        if (Y_boot$type=="distribution"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w,,drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }

        if (Y_boot$type=="sphere"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w,,drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }

        if (Y_boot$type=="sphere2"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w,,drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }

        if (Y_boot$type=="covariance"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[,,w, drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }
        
        F_SPLIT <- NULL
        decoupe <- 0
        
        if (is.element("factor",split.spaces)==TRUE){
          
          if( ERT==FALSE){
            feuille_split_Factor <- var_split(Factor_courant,Y_courant,timeScale)
          }
          
          else{feuille_split_Factor <- ERvar_split(Factor_courant,Y_courant,timeScale,ntry = ntry)}
          
          if (feuille_split_Factor$Pure==FALSE){
            F_SPLIT <- rbind(F_SPLIT,c("Factor",feuille_split_Factor$impurete))
            decoupe <- decoupe +1
          }
        }
        
        if (is.element("curve",split.spaces)==TRUE){
          
          if( ERT==FALSE){
            feuille_split_Curve <- var_split(Curve_courant,Y_courant,timeScale)
          }
          
          else{feuille_split_Curve <- ERvar_split(Curve_courant,Y_courant,timeScale, ntry=ntry)}
          
          if (feuille_split_Curve$Pure==FALSE){
            F_SPLIT <- rbind(F_SPLIT,c("Curve",feuille_split_Curve$impurete))
            decoupe <- decoupe +1
          }
        }
        
        if (is.element("scalar",split.spaces)==TRUE){
          
          if( ERT==FALSE){
            feuille_split_Scalar <- var_split(Scalar_courant,Y_courant,timeScale)
          }
          
          else{feuille_split_Scalar <- ERvar_split(Scalar_courant,Y_courant,timeScale, ntry=ntry)}
          
          if (feuille_split_Scalar$Pure==FALSE){
            F_SPLIT <- rbind(F_SPLIT,c("Scalar",feuille_split_Scalar$impurete))
            decoupe <- decoupe +1
          }
        }
        
        if (is.element("shape",split.spaces)==TRUE){
          
          feuille_split_Shape <- ERvar_split(Shape_courant,Y_courant,timeScale, ntry=ntry)
          
          
          if (feuille_split_Shape$Pure==FALSE){
            F_SPLIT <- rbind(F_SPLIT,c("Shape",feuille_split_Shape$impurete))
            decoupe <- decoupe +1
          }
          
        }
        
        
        if (is.element("image",split.spaces)==TRUE){
          
          feuille_split_Image <- ERvar_split(Image_courant,Y_courant,timeScale, ntry=ntry)
          
          if (feuille_split_Image$Pure==FALSE){
            F_SPLIT <- rbind(F_SPLIT,c("Image",feuille_split_Image$impurete))
            decoupe <- decoupe +1
          }
          
        }
        
        
        if (decoupe>0){
          
          TYPE <- F_SPLIT[which.min(F_SPLIT[,2]),1]
          X <- get(TYPE)
          X_boot <- get(paste(TYPE,"_boot",sep=""))
          
          feuille_split <- get(paste("feuille_split_",TYPE, sep=""))
          
          vsplit_space <- get(paste("tirage",TYPE, sep=""))[feuille_split$variable]
          
          #if (imp_apres_split<imp_avant_split){
          
          gauche_id <- unique(Y_boot$id[w])[which(feuille_split$split==1)]
          droit_id <- unique(Y_boot$id[w])[which(feuille_split$split==2)]
          
          if (Y$type=="surv"){
            imp_nodes[[2*unique(id_feuille)[i]]] <- Inf
            imp_nodes[[2*unique(id_feuille)[i]+1]] <- Inf
          }
          else {
            imp_nodes[[2*unique(id_feuille)[i]]] <- feuille_split$impur_list[[1]]
            imp_nodes[[2*unique(id_feuille)[i]+1]] <- feuille_split$impur_list[[2]]
          }
          
          
          V_split <- rbind(V_split,c(TYPE,unique(id_feuille)[i],vsplit_space))
          
          wY_gauche <- NULL
          wY_droit <- NULL
          w_gauche <- NULL
          w_droit <- NULL
          
          
          for (k in 1:length(gauche_id)){
            w_gauche <- c(w_gauche, which(X_boot$id==gauche_id[k]))
            wY_gauche <- c(wY_gauche, which(Y_boot$id==gauche_id[k]))
          }
          
          for (k in 1:length(droit_id)){
            w_droit <- c(w_droit, which(X_boot$id==droit_id[k]))
            wY_droit <- c(wY_droit, which(Y_boot$id==droit_id[k]))
          }
          
          
          id_feuille_prime[wY_gauche] <- 2*(unique(id_feuille)[i])
          id_feuille_prime[wY_droit] <- 2*(unique(id_feuille)[i])+1
          
          #print(paste("Split on the variable", vsplit_space, "on the space of ", paste(TYPE,"s",sep="")))
          
          if (X$type=="curve"){
            trajG <- as.data.frame(cbind(X_boot$id[w_gauche], X_boot$time[w_gauche], X_boot$X[w_gauche,vsplit_space]))
            trajD <- as.data.frame(cbind(X_boot$id[w_droit], X_boot$time[w_droit], X_boot$X[w_droit,vsplit_space]))
            meanFg <- as.matrix(kmlShape::meanFrechet(trajG))
            meanFd <- as.matrix(kmlShape::meanFrechet(trajD))
          }
          
          if (X$type=="shape"){
            meanFg <- mshape(X_boot$X[,,w_gauche,vsplit_space])
            meanFd <- mshape(X_boot$X[,,w_droit,vsplit_space])
          }
          
          if (X$type=="image"){
            meanFg <- apply(X_boot$X[w_gauche,,vsplit_space, drop=FALSE],2,"mean")
            meanFd <- apply(X_boot$X[w_droit,,vsplit_space, drop=FALSE],2,"mean")
          }
          
          if (X$type=="factor"){
            meanFg <- unique(X_boot$X[w_gauche, vsplit_space])
            meanFd <- unique(X_boot$X[w_droit,vsplit_space])
          }
          
          if (X$type=="scalar"){
            meanFg <- mean(X_boot$X[w_gauche,vsplit_space])
            meanFd <- mean(X_boot$X[w_droit,vsplit_space])
          }
          
          
          hist_nodes[[2*(unique(id_feuille)[i])]] <- meanFg
          hist_nodes[[2*(unique(id_feuille)[i])+1]] <- meanFd
          count_split <- count_split+1
          
          feuilles_courantes <- unique(id_feuille_prime)
        }
        
        
      }
    }
    
    id_feuille <- id_feuille_prime
    
    if (count_split ==0 ){
      
      V_split <- data.frame(V_split)
      names(V_split) <- c("type","num_noeud", "var_split")
      for (q in unique(id_feuille)){
        w <- which(id_feuille == q)
        if (Y$type=="curve"){
          Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y_boot$id[w], Y_boot$time[w], Y_boot$Y[w]))
        }
        if(Y$type=="scalar"){
          Y_pred[[q]]<- mean(Y_boot$Y[w])
        }
        if(Y$type=="factor"){
          Table <- which.max(table(Y_boot$Y[w]))
          Y_pred[[q]] <-  as.factor(attributes(Table)$names)
        }
        
        if (Y$type=="shape"){
          Y_pred[[q]] <-  mshape(Y_boot$Y[,,w, drop=FALSE])
        }
        
        if (Y$type=="image"){
          Y_pred[[q]] <- apply(Y_boot$Y[w,,drop=FALSE], 2, "mean")
        }
        
        if (Y$type=="surv"){
          donnees = survfit(Surv(Y_boot$time[w], Y_boot$Y[w])~1)
          Y_pred[[q]] <- data.frame(times=donnees$time, traj=donnees$surv)
          Y_pred_surv[[q]] <- data.frame(time=Y_boot$time[w], traj = Y_boot$Y[w])
        }
        
        if (Y$type=="distribution"){
          Y_pred[[q]] <- colMeans(Y_boot$Y[w,,drop=FALSE])
        }

        if (Y$type=="sphere"){
          S <- list()
          nRow <- nrow(Y_boot$Y[w,,drop=FALSE])
          for (h in 1:nRow){
            S[[h]] = Y_boot$Y[w,,drop=FALSE][h,] # project onto Sphere
          }
          S = riemfactory(S, name="sphere")
          Y_pred[[q]] <- as.vector(rbase.mean(S)$x)
        }

        if (Y$type=="sphere2"){
          Y_pred[[q]] <- IntrinsicMean(Y_boot$Y[w,,drop=FALSE])
        }

        if (Y$type=="covariance"){
          Y_pred[[q]] <- estcov(Y_boot$Y[,,w, drop=FALSE], method= method, alpha=alpha)$mean
        }         
        
      }
      if (Y$type=="factor"){
        Ylevels <- unique(Y_boot$Y)
        return(list(feuilles = id_feuille, idY=Y_boot$id,Ytype=Y_boot$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, time = time, Y=Y, boot=boot, Ylevels=Ylevels))
      }
      return(list(feuilles = id_feuille, idY=Y_boot$id,Ytype=Y_boot$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, time = time, Y=Y, boot=boot, Y_pred_surv=Y_pred_surv))
    }
  }
  
  
  V_split <- data.frame(V_split)
  names(V_split) <- c("type","num_noeud", "var_split")
  for (q in unique(id_feuille)){
    w <- which(id_feuille == q)
    if (Y$type=="curve"){
      Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y_boot$id[w], Y_boot$time[w], Y_boot$Y[w]))
    }
    
    if (Y$type=="image"){
      Y_pred[[q]] <- apply(Y_boot$Y[w,,drop=FALSE], 2, "mean")
    }
    
    if(Y$type=="scalar"){
      Y_pred[[q]]<- mean(Y_boot$Y[w])
    }
    
    if(Y$type=="factor"){
      Table <- which.max(table(Y_boot$Y[w]))
      Y_pred[[q]] <-  as.factor(attributes(Table)$names)
    }
    
    if (Y$type=="shape"){
      Y_pred[[q]] <- mshape(Y_boot$Y[,,w, drop=FALSE])
    }
    
    if (Y$type=="surv"){
      donnees = survfit(Surv(Y_boot$time[w], Y_boot$Y[w])~1)
      Y_pred[[q]] <- data.frame(times=donnees$time, traj=donnees$surv)
      Y_pred_surv[[q]] <- data.frame(time=Y_boot$time[w], traj = Y_boot$Y[w])
    }
    
    if (Y$type=="distribution"){
      Y_pred[[q]] <- colMeans(Y_boot$Y[w,,drop=FALSE])
    }  

    if (Y$type=="sphere"){
      S <- list()
      nRow <- nrow(Y_boot$Y[w,,drop=FALSE])
      for (h in 1:nRow){
        S[[h]] = Y_boot$Y[w,,drop=FALSE][h,] # project onto Sphere
      }
      S = riemfactory(S, name="sphere")
      Y_pred[[q]] <- as.vector(rbase.mean(S)$x)
    }

    if (Y$type=="sphere2"){
      Y_pred[[q]] <- IntrinsicMean(Y_boot$Y[w,,drop=FALSE])
    }

    if (Y$type=="covariance"){
      Y_pred[[q]] <- estcov(Y_boot$Y[,,w, drop=FALSE], method= method, alpha=alpha)$mean
    }   
    
  }
  if (Y$type=="factor"){
    Ylevels <- unique(Y_boot$Y)
    return(list(feuilles = id_feuille, idY=Y_boot$id,Ytype=Y_boot$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, time = time, Y=Y,boot=boot, Ylevels=Ylevels))
  }
  return(list(feuilles = id_feuille, idY=Y_boot$id, Ytype=Y_boot$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred= Y_pred, time=time, Y=Y, boot=boot, Y_pred_surv=Y_pred_surv))
}


#' Predict with Frechet random forests
#'
#' @param object : Frechet random forest
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Shape [list]:
#' @param Image [list]:
#' @param aligned.shape [logical]:
#' @param timeScale [numeric]:
#' @param d_out [numeric]:
#' @param ... : optional parameters to be passed to the low level function
#'
#' @import kmlShape
#' @import stringr
#' @import Evomorph
#' @import geomorph
#' @import survival
#'
#' @export
#'
predict.FrechForest <- function(object, Curve=NULL,Scalar=NULL,Factor=NULL,Shape=NULL, Image=NULL,aligned.shape=FALSE, timeScale=0.1, d_out=0.1,...){
  # La première etape est de toujours lire les predicteurs ::
  
  if (is.null(Curve)==FALSE){
    Curve <- list(type="curve",X=Curve$X,id=Curve$id,time=Curve$time)
  }
  if (is.null(Scalar)==FALSE){
    Scalar <- list(type="scalar",X=Scalar$X,id=Scalar$id)
  }
  if (is.null(Factor)==FALSE){
    Factor <- list(type="factor",X=Factor$X,id=Factor$id)
  }
  if (is.null(Shape)==FALSE){
    Shape <- list(type="shape",X=Shape$X,id=Shape$id)
  }
  if (is.null(Image)==FALSE){
    Image <- list(type="image",X=Image$X,id=Image$id)
  }
  
  ## Puis on prend les predicteurs:
  
  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs
  # On va les lires en mettant la maj sur les differents elements qui le constituent :
  
  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }
  
  if (is.element("shape",inputs)==TRUE & aligned.shape==FALSE){
    for (k in 1:dim(Shape$X)[length(dim(Shape$X))]){
      Shape$X[,,,k] <- gpagen(Shape$X[,,,k],print.progress = FALSE)$coords
    }
    aligned.shape=TRUE
  }
  
  Id.pred <- unique(get(Inputs[1])$id)
  pred.feuille <- matrix(0, ncol(object$rf), length(Id.pred))
  
  if (object$type=="factor"){
    pred.feuille <- as.data.frame(matrix(0, ncol(object$rf), length(Id.pred)))
  }
  
  for (t in 1:ncol(object$rf)){
    pred.feuille[t,] <- pred.FT(object$rf[,t], Curve = Curve,Scalar = Scalar,Factor=Factor,Shape=Shape,Image=Image, timeScale, aligned.shape = aligned.shape)
  }
  
  if (object$type=="scalar"){
    pred <- apply(pred.feuille, 2, "mean")
    return(pred)
  }
  
  if (object$type=="factor"){
    pred.all <- apply(pred.feuille, 2, "table")
    val <- factor(rep(NA, length(pred.all)), levels=object$levels)
    proba <- rep(NA, length(pred.all))
    for (k in 1:length(pred.all)){
      val[k] <- factor(attributes(which.max(pred.all[[k]])))
      proba[k] <- max(pred.all[[k]])/ncol(object$rf)
    }
    prediction <- data.frame(pred=val, prob=proba)
    return(prediction)
  }
  
  if (object$type=="shape"){
    pred <- array(0, dim=c(object$size[1], object$size[2],length(Id.pred)))
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- array(0,dim=c(object$size[1],object$size[2],ncol(object$rf)))
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant[,,k] <- object$rf[,k]$Y_pred[[pred.feuille[k,l]]]
      }
      Ms <- mshape(pred_courant)
      M <- matrix(0,dim(Ms)[1], dim(Ms)[2])
      M[,1] <- Ms[,1]
      M[,2] <- Ms[,2]
      pred[,,l] <- M
    }
  }
  
  if (object$type=="image"){
    pred <- matrix(0, length(Id.pred), object$size[2])
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- matrix(0,ncol(object$rf),object$size[2])
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant[,,k] <- object$rf[,k]$Y_pred[[pred.feuille[k,l]]]
      }
      pred[,,l] <-  apply(pred_courant, 2, "mean")
    }
  }
  
  if (object$type=="distribution"){
    pred <- matrix(0, length(Id.pred), object$size[2])
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- matrix(0,ncol(object$rf),object$size[2])
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant[k,] <- object$rf[,k]$Y_pred[[pred.feuille[k,l]]]
      }
      pred[l,] <-  colMeans(pred_courant)
    }
  }

  if (object$type=="sphere"){
    pred <- matrix(0, length(Id.pred), object$size[2])
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- matrix(0,ncol(object$rf),object$size[2])
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant[k,] <- object$rf[,k]$Y_pred[[pred.feuille[k,l]]]
      }
      S <- list()
      nRow <- nrow(pred_courant)
      for (h in 1:nRow){
        S[[h]] = pred_courant[h,] # project onto Sphere
      }
      S = riemfactory(S, name="sphere")  
      pred[l,] <- as.vector(rbase.mean(S)$x)
    }
  }

  if (object$type=="sphere2"){
    pred <- matrix(0, length(Id.pred), object$size[2])
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- matrix(0,ncol(object$rf),object$size[2])
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant[k,] <- object$rf[,k]$Y_pred[[pred.feuille[k,l]]]
      }
      pred[l,] <-  IntrinsicMean(pred_courant)
    }
  }

  if (object$type=="covariance"){
    pred <- array(0, dim=c(object$size[1], object$size[2],length(Id.pred)))
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- array(0,dim=c(object$size[1],object$size[2],ncol(object$rf)))
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant[,,k] <- object$rf[,k]$Y_pred[[pred.feuille[k,l]]]
      }
      pred[,,l] <-  estcov(pred_courant, method= method, alpha=alpha)$mean
    }
  }
  
  if (object$type=="curve"){
    pred <- NULL
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- NULL
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant <- rbind(pred_courant, cbind(rep(k,dim(object$rf[,k]$Y_pred[[pred.feuille[k,l]]])[1]),object$rf[,k]$Y_pred[[pred.feuille[k,l]]]))
      }
      predouille <- kmlShape::meanFrechet(pred_courant, timeScale = timeScale)
      predouille <- cbind(predouille, rep(Id.pred[l],dim(predouille)[1]))
      pred <- rbind(pred, predouille)
    }
    names(pred) <- c("times", "traj", "ID")
  }
  
  if (object$type=="surv"){
    pred <- NULL
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- NULL
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant <- rbind(pred_courant, object$rf[,k]$Y_pred_surv[[pred.feuille[k,l]]])
      }
      predouille <- survfit(Surv(pred_courant$time[w], pred_courant$traj[w])~1)
      P = data.frame(times=predouille$time, traj=predouille$surv, id= rep(l,length(predouille$time)))
      pred <- rbind(pred, P)
    }
    #names(pred) <- c("times", "traj", "ID")
  }
  return(list(pred=pred,pred_feuille=pred.feuille))
}

#' OOB tree
#'
#' @param tree
#' @param Curve
#' @param Scalar
#' @param Factor
#' @param Shape
#' @param Image
#' @param Y
#' @param timeScale
#' @param d_out
#'
#' @import kmlShape
#' @import Evomorph
#' @import stringr
#'
#' @keywords internal
OOB.tree <- function(tree, Curve=NULL, Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL ,Y, timeScale=0.1, d_out=0.1){
  
  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs
  
  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }
  
  BOOT <- tree$boot
  OOB <- setdiff(unique(Y$id), BOOT)
  xerror <- rep(NA,length(OOB))
  Scalar_courant <- NULL
  Factor_courant <- NULL
  Curve_courant <- NULL
  Image_courant <- NULL
  Shape_courant <- NULL
  if (Y$type=="curve" || Y$type=="surv"){
    for (i in OOB){
      id_wY <- which(Y$id== i)
      if (is.element("curve",inputs)==TRUE) {
        id_wXCurve <- which(Curve$id==i)
        Curve_courant <- list(type="curve",X=Curve$X[id_wXCurve,,drop=FALSE], id=Curve$id[id_wXCurve],time=Curve$time[id_wXCurve])
      }
      
      if (is.element("shape",inputs)==TRUE){
        id_wXShape <- which(Shape$id==i)
        Shape_courant <- list(type="shape",X=Shape$X[,,id_wXShape,,drop=FALSE], id=Shape$id[id_wXShape])
      }
      
      if (is.element("image",inputs)==TRUE){
        id_wXImage <- which(Image$id==i)
        Image_courant <- list(type="image",X=Image$X[id_wXImage,,,drop=FALSE], id=Image$id[id_wXImage])
      }
      
      if (is.element("factor",inputs)==TRUE){
        id_wXFactor <- which(Factor$id==i)
        Factor_courant <- list(type="factor",X=Factor$X[id_wXFactor,,drop=FALSE], id=Factor$id[id_wXFactor])
      }
      
      if (is.element("scalar",inputs)==TRUE){
        id_wXScalar <- which(Scalar$id==i)
        Scalar_courant <- list(type="scalar",X=Scalar$X[id_wXScalar,,drop=FALSE], id=Scalar$id[id_wXScalar])
      }
      pred_courant <- pred.FT(tree, Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape = Shape_courant,Image=Image_courant, timeScale=timeScale, aligned.shape = TRUE)
      #chancla <- DouglasPeuckerNbPoints(tree$Y_Curves[[pred_courant]]$times, tree$Y_Curves[[pred_courant]]$traj, nbPoints = length(stats::na.omit(Y[id_w])))
      xerror[which(OOB==i)] <- kmlShape::distFrechet(tree$Y_pred[[pred_courant]]$times, tree$Y_pred[[pred_courant]]$traj, Y$time[id_wY], Y$Y[id_wY], timeScale = d_out)^2
    }
  }
  else {
    w_XCurve <- NULL
    w_XScalar <- NULL
    w_XFactor <- NULL
    w_XShape <- NULL
    w_XImage <- NULL
    w_y <- NULL
    for (i in OOB){
      
      if (is.element("curve",inputs)==TRUE) w_XCurve <- c(w_XCurve, which(Curve$id==i))
      if (is.element("scalar",inputs)==TRUE) w_XScalar <- c(w_XScalar, which(Scalar$id==i))
      if (is.element("factor",inputs)==TRUE) w_XFactor <- c(w_XFactor, which(Factor$id==i))
      if (is.element("shape",inputs)==TRUE) w_XShape <- c(w_XShape, which(Shape$id==i))
      if (is.element("image",inputs)==TRUE) w_XImage <- c(w_XImage, which(Image$id==i))
      
      w_y <- c(w_y, which(Y$id==i))
    }
    
    if (is.element("curve",inputs)==TRUE) Curve_courant <- list(type="curve",X=Curve$X[w_XCurve,,drop=FALSE], id=Curve$id[w_XCurve],time=Curve$time[w_XCurve])
    if (is.element("scalar",inputs)==TRUE) Scalar_courant  <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
    if (is.element("factor",inputs)==TRUE) Factor_courant  <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
    if (is.element("shape",inputs)==TRUE) Shape_courant  <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
    if (is.element("image",inputs)==TRUE) Image_courant  <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
    
    pred <- pred.FT(tree,Curve=Curve_courant,Scalar = Scalar_courant,Factor=Factor_courant, Shape=Shape_courant, Image = Image_courant, aligned.shape = TRUE)
    
    if (Y$type=="scalar"){xerror <- (Y$Y[w_y]-pred)^2}
    if (Y$type=="factor"){xerror <- 1*(pred!=Y$Y[w_y])}
    
    if (Y$type=="shape"){
      xerror <- rep(NA,length(pred))
      for (l in 1:length(pred)){
        xerror[l] <- ShapeDist(Y$Y[,,w_y[l], drop=FALSE],tree$Y_pred[[pred[l]]])^2
      }
    }
    
    if (Y$type=="image"){
      xerror <- rep(NA,length(pred))
      for (l in 1:length(pred)){
        xerror[l] <- mean((Y$Y[w_y[l],]-tree$Y_pred[[pred[l]]])^2)
      }
    }
    
    if (Y$type=="distribution"){
      xerror <- rep(NA,length(pred))
      for (l in 1:length(pred)){
        xerror[l] <- pracma::trapz(qSup, (Y$Y[w_y[l],] - tree$Y_pred[[pred[l]]])^2)
      }
    }

    if (Y$type=="sphere"){
      xerror <- rep(NA,length(pred))
      for (l in 1:length(pred)){
        tmp = crossprod(Y$Y[w_y[l],],tree$Y_pred[[pred[l]]])[1,1]
        tmp = sign(tmp) * min(1, abs(tmp))
        xerror[l] <- (acos(tmp))^2
      }
    }

    if (Y$type=="sphere2"){
      xerror <- rep(NA,length(pred))
      for (l in 1:length(pred)){
        xerror[l] <- distCosine(Y$Y[w_y[l],],tree$Y_pred[[pred[l]]],r=1)^2
      }
    }

    if (Y$type=="covariance"){
      xerror <- rep(NA,length(pred))
      for (l in 1:length(pred)){
        xerror[l] <- distcov(Y$Y[,,w_y[l]], tree$Y_pred[[pred[l]]], method=method, alpha=alpha)^2
      }
    }
    
  }
  return(mean(xerror))
}

#' Title
#'
#' @param time.init
#' @param traj.init
#' @param time.new
#'
#'
#' @keywords internal
Curve.reduc.times <- function(time.init , traj.init, time.new){
  new.Curve <- matrix(NA,length(time.new),2)
  for (j in 1:length(time.new)){
    w.time <- which.min(abs(time.new[j]-time.init))
    if (round(time.init[w.time]-time.new[j],5)==0){
      new.Curve[j,] <- c(time.new[j], traj.init[w.time])
    }
    else {
      t_g <- (time.new[j]>time.init[w.time])*(time.init[w.time]) + (time.new[j]<time.init[w.time])*(time.init[w.time-1])
      t_d <- (time.new[j]<time.init[w.time])*(time.init[w.time]) + (time.new[j]>time.init[w.time])*(time.init[w.time+1])
      Y_g <- (time.new[j]>time.init[w.time])*(traj.init[w.time]) + (time.new[j]<time.init[w.time])*(traj.init[w.time-1])
      Y_d <- (time.new[j]<time.init[w.time])*(traj.init[w.time]) + (time.new[j]>time.init[w.time])*(traj.init[w.time+1])
      d1 <- time.new[j]-t_g
      d2 <- t_d - time.new[j]
      new.Curve[j,] <- c(time.new[j], (1 - (d1/(d1+d2)))*Y_g + (1 - (d2/(d1+d2)))*Y_d)
    }
  }
  return(new.Curve)
}

#' OOB for random Forest
#'
#' @param rf
#' @param Curve
#' @param Scalar
#' @param Factor
#' @param Shape
#' @param Image
#' @param Y
#' @param timeScale
#' @param d_out
#'
#' @import stringr
#' @import kmlShape
#' @import Evomorph
#' @import geomorph
#'
#' @keywords internal
OOB.rfshape <- function(rf, Curve=NULL, Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL, Y, timeScale=0.1, d_out=0.1){
  
  ### Pour optimiser le code il faudra virer cette ligne et ne le calculer qu'une seule fois !
  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs
  
  
  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }
  
  
  err <- rep(NA,length(unique(Y$id)))
  
  Curve_courant <- NULL
  Scalar_courant <- NULL
  Factor_courant <- NULL
  Shape_courant <- NULL
  Image_courant <- NULL
  
  
  if (Y$type=="surv"){
    oob.pred <- list()
    #errdp <- rep(NA,length(unique(id)))
    
    for (i in 1:length(unique(Y$id))){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      courbe2 = Y$Y[w_y]
      pred_courant <- rep(0, length(w_y))
      n_ind = 0
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){
          n_ind= n_ind+1
          
          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }
          
          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }
          
          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }
          
          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }
          
          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }
          
          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale, aligned.shape = TRUE)
          courbe <- rf$rf[,t]$Y_pred[[pred]] ## Il faut les mettre aux mêmes temps que Y$time[w_y]
          
          for (j in 1:length(w_y)){
            courbe2[j] = courbe[which.min(abs(Y$time[w_y][j]-courbe[,1])),2]
          }
          pred_courant <- pred_courant + courbe2
        }
      }
      oob.pred[[i]] <-  data.frame(times=Y$time[w_y], traj=pred_courant/n_ind)
      err[i] <- mean((Y$Y[w_y]-pred_courant/n_ind)^2)
    }
    return(list(err=err,oob.pred=oob.pred))
  }
  
  if (Y$type=="curve"){
    oob.pred <- list()
    #errdp <- rep(NA,length(unique(id)))
    
    for (i in 1:length(unique(Y$id))){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){
          
          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }
          
          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }
          
          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }
          
          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }
          
          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }
          
          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale, aligned.shape = TRUE)
          courbe <- rf$rf[,t]$Y_pred[[pred]]
          pred_courant <- rbind(cbind(rep(t,dim(courbe)[1]),courbe),pred_courant)
        }
      }
      mean_pred <- meanFrechet(pred_courant, timeScale = d_out)
      dp <- as.data.frame(Curve.reduc.times(mean_pred$times, mean_pred$traj, Y$time[w_y]))
      names(dp) <- c("x","y")
      oob.pred[[i]] <- dp
      err[i] <- distFrechet(dp$x, dp$y, Y$time[w_y], Y$Y[w_y], timeScale = d_out)^2
    }
    return(list(err=err,oob.pred=oob.pred))
  }
  
  if (Y$type=="scalar"){
    oob.pred <- rep(NA, length(unique(Y$id)))
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- Y$id[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){
          
          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }
          
          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }
          
          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }
          
          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }
          
          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }
          
          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale, aligned.shape = TRUE)
          pred_courant <- c(pred_courant, pred)
        }
      }
      oob.pred[i] <- mean(pred_courant)
      err[i] <- (oob.pred[i]-Y$Y[w_y])^2
    }
  }
  
  if (Y$type=="factor"){
    oob.pred <- factor(rep(NA, length(unique(Y$id))), levels=rf$levels)
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- Y$id[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- factor(rep(NA, length(unique(Y$id))), levels=rf$levels)
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){
          
          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }
          
          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }
          
          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }
          
          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }
          
          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }
          
          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale, aligned.shape = TRUE)
          pred_courant[t] <- pred
        }
      }
      pred_courant <- na.omit(pred_courant)
      oob.pred[i] <- as.factor(attributes(which.max(table(pred_courant))))
    }
    err <- 1*(oob.pred!=Y$Y)
  }
  
  if (Y$type=="shape"){
    oob.pred <- array(0,dim=dim(Y$Y))
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- array(0, dim=c(dim(Y$Y)[1],dim(Y$Y)[2],length(rf$rf)))
      selection <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){
          
          selection <- c(selection, t)
          
          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }
          
          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }
          
          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }
          
          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }
          
          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }
          
          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale, aligned.shape = TRUE)
          pred_courant[,,t] <- rf$rf[,t]$Y_pred[[pred]]
        }
      }
      pred_courant <- pred_courant[,,selection]
      mean_pred <- mshape(pred_courant)
      err[i] <- ShapeDist(Y$Y[,,w_y, drop=FALSE],mean_pred)^2
      M <- matrix(0,dim(mean_pred)[1], dim(mean_pred)[2])
      M[,1] <- mean_pred[,1]
      M[,2] <- mean_pred[,2]
      
      oob.pred[,,i] <- M
    }
    return(list(err=err,oob.pred=oob.pred))
  }
  
  if (Y$type=="image"){
    oob.pred <- array(0,dim=dim(Y$Y))
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- array(0, dim=c(length(rf$rf),ncol(Y$Y)))
      selection <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){
          selection <- c(selection, t)
          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }
          
          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }
          
          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }
          
          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }
          
          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }
          
          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale, aligned.shape = TRUE)
          pred_courant[t,] <- rf$rf[,t]$Y_pred[[pred]]
        }
      }
      oob.pred[i,] <-  apply(pred_courant,2,"mean")
      err[i] <- mean((oob.pred[i,]-Y$Y[w_y,])^2)
    }
    return(list(err=err,oob.pred=oob.pred))
  }
  
  if (Y$type=="distribution"){
    oob.pred <- array(NA,dim=dim(Y$Y))
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- NULL
      #selection <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        #print(oob)
        if (is.element(indiv, oob)== TRUE){
          #selection <- c(selection, t)
          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }
          
          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }
          
          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }
          
          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }
          
          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }
          
          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale, aligned.shape = TRUE)
          pred_courant <- rbind(pred_courant, rf$rf[,t]$Y_pred[[pred]]) 
          
        }
      }
      #print(pred_courant)
      if (!is.null(pred_courant)){
        mean_pred <- colMeans(qin = pred_courant)
        oob.pred[i,] <- mean_pred
        err[i] <- pracma::trapz(qSup, (mean_pred - Y$Y[w_y,])^2)
      }
    }
    return(list(err=err,oob.pred=oob.pred))
  }

  if (Y$type=="sphere"){
    oob.pred <- array(NA,dim=dim(Y$Y))
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- NULL
      #selection <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        #print(oob)
        if (is.element(indiv, oob)== TRUE){
          #selection <- c(selection, t)
          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }
          
          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }
          
          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }
          
          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }
          
          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }
          
          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale, aligned.shape = TRUE)
          pred_courant <- rbind(pred_courant, rf$rf[,t]$Y_pred[[pred]]) 
          
        }
      }
      #print(pred_courant)
      if (!is.null(pred_courant)){
        S <- list()
        nRow <- nrow(pred_courant)
        for (h in 1:nRow){
          S[[h]] = pred_courant[h,] # project onto Sphere
        }
        S = riemfactory(S, name="sphere")  
        mean_pred <- as.vector(rbase.mean(S)$x)
        oob.pred[i,] <- mean_pred 
        tmp = crossprod(mean_pred,Y$Y[w_y,])[1,1]
        tmp = sign(tmp) * min(1, abs(tmp))
        err[i] <- (acos(tmp))^2
      }
    }
    return(list(err=err,oob.pred=oob.pred))
  }

  if (Y$type=="sphere2"){
    oob.pred <- array(NA,dim=dim(Y$Y))
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- NULL
      #selection <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        #print(oob)
        if (is.element(indiv, oob)== TRUE){
          #selection <- c(selection, t)
          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }
          
          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }
          
          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }
          
          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }
          
          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }
          
          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale, aligned.shape = TRUE)
          pred_courant <- rbind(pred_courant, rf$rf[,t]$Y_pred[[pred]]) 
          
        }
      }
      #print(pred_courant)
      if (!is.null(pred_courant)){
        mean_pred <- IntrinsicMean(qin = pred_courant)
        oob.pred[i,] <- mean_pred
        err[i] <- distCosine(mean_pred,Y$Y[w_y,],r=1)^2
      }
    }
    return(list(err=err,oob.pred=oob.pred))
  }

  if (Y$type=="covariance"){
    oob.pred <- array(NA,dim=dim(Y$Y))
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- array(0, dim=c(dim(Y$Y)[1],dim(Y$Y)[2],length(rf$rf)))
      selection <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){
          
          selection <- c(selection, t)
          
          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve])
          }
          
          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }
          
          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }
          
          if (is.element("shape",inputs)==TRUE){
            w_XShape <- which(Shape$id== indiv)
            Shape_courant <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
          }
          
          if (is.element("image",inputs)==TRUE){
            w_XImage <- which(Image$id== indiv)
            Image_courant <- list(type="image", X=Image$X[w_XImage,,, drop=FALSE], id=Image$id[w_XImage])
          }
          
          pred <- pred.FT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape=Shape_courant,Image=Image_courant, timeScale = timeScale, aligned.shape = TRUE)
          pred_courant[,,t] <- rf$rf[,t]$Y_pred[[pred]]
        }
      }
      if (!is.null(selection)){
        pred_courant <- pred_courant[,,selection,drop=FALSE]
        mean_pred <- estcov(pred_courant, method= method, alpha=alpha)$mean
        err[i] <- distcov(Y$Y[,,w_y], mean_pred, method=method, alpha=alpha)^2
        oob.pred[,,i] <- mean_pred
      }
    }
    return(list(err=err,oob.pred=oob.pred))
  }
  return(list(err=err,oob.pred=oob.pred))
}

#' Title
#'
#' @param Courbes
#' @param id
#'
#' @keywords internal
permutation_courbes <- function(Courbes,id){
  perm <- sample(unique(id), length(unique(id))) #### on change les identifiants ::
  new <- NULL
  for (i in perm){
    w <- which(id==i)
    new <- c(new,Courbes[w])
  }
  return(new)
}

#' Title
#'
#' @param Shapes
#' @param id
#'
#' @keywords internal
permutation_shapes <- function(Shapes, id){
  perm <- sample(id,length(id))
  new <- array(0,dim=dim(Shapes)[1:3])
  for (i in 1:length(id)){
    new[,,i] <- Shapes[,,which(id==perm[i])]
  }
  return(new)
}

#' Frechet Random Forest
#'
#' This function builds Frechet random Forest introduced by Capitaine et.al, this includes the 
#, OOB errors and variable importance computations.
#'
#'
#' @param Curve [list]: A list that contains the different input curves. It must contain the following elements (no choice): \code{X} the matrix of the different curves, each column code for a different curve variable; \code{id} is the vector of the identifiers for the different trajectories contained in \code{X}; \code{time} is the vector of the measurement times associated with the trajectories contained in \code{X}.
#' @param Scalar [list]: A list that contains the different input scalars. It must contain the following elements (no choice):  \code{X} the matrix of the scalars, each column code for a different variable; \code{id} is the vector of the identifiers for each individual.
#' @param Factor [list]: A list that contains the different input factors. It must contain the following elements (no choice):  \code{X} the matrix of the factors, each column code for a different variable; \code{id} is the vector of the identifiers for each individual.
#' @param Shape [list]: A list that contains the different input shapes. It must contain the following elements (no choice):  \code{X} the array of the shapes of dimension \code{n}x2x\code{l}x\code{p} where \code{n} is the number of points for composing each shape, \code{l} is the number of shapes and \code{p} is the number of shapes variables, \code{id} is the vector of the identifiers for each individual.
#' @param Image [list]: A list that contains the different input images. It must contain the following elements (no choice):  \code{X} the array of the images of dimension \code{n}x\code{m}x\code{l}x\code{p} where \code{n}*\code{m} is the size of each image, \code{l} is the number of images and \code{p} is the number of shapes variables; \code{id} is the vector of the identifiers for each individual.
#' @param Y [list]: A list that contains the output, It must contain the following elements (no choice): \code{type} defines the nature of the output, can be "\code{curve}", "\code{sclalar}", "\code{factor}", "\code{shape}", "\code{image}"; \code{Y} is the output variable; \code{id} is the vector of the identifiers for each individuals, they should be the same as the identifiers of the inputs.
#' @param mtry [numeric]: Number of variables randomly sampled as candidates at each split. The default value \code{p/3}
#' @param ntree [numeric]: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
#' @param ncores [numeric]: Number of cores used to build Frechet randomized trees in parallel, defaulting to number of cores of the computer minus 1.
#' @param ERT [logical]: If \code{TRUE} uses Extremly Randomized Frechet Trees to build the Frechet forest.
#' @param ntry [numeric]: Only with \code{ERT=TRUE}, allows to manage with randomness of the trees.
#' @param timeScale [numeric]: Allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping. Only used when there are trajectories either in input or output.
#' @param imp [logical]: TRUE to compute the variables importance FALSE otherwise (default \code{imp=}TRUE)
#' @param d_out [string]: "euc" or "frec".
#' @param ... : optional parameters to be passed to the low level function
#'
#' @import stringr
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import pbapply
#'
#' @return A Frechet random forest which is a list of the following elements: \itemize{
#' \item \code{rf:} a list of the \code{ntree} randomized Frechet trees that compose the forest.
#' \item \code{xerror :} a vector containing the OOB prediction error of each randomized Frechet tree composing the forest.
#' \item \code{OOB.err: } a vector containing the OOB prediction error of each individual in the learning sample.
#' \item \code{OOB.pred: } a list of the OOB prediction for each individual in the learning set.
#' \item \code{Importance: } A vector containing the variables importance.
#' \item \code{varex: } “pseudo R-squared”: Percentage of variance explained.
#' }
#' @export
#'
FrechForest <- function(Curve=NULL,Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL ,Y, mtry=NULL, ntree=100,ncores=NULL,ERT=FALSE, timeScale=0.1,ntry=3, OOB=TRUE, imp=TRUE, d_out=0.1, ...){
  
  
  ### On va regarder les differentes entrees:
  if (is.null(Curve)==FALSE){
    Curve <- list(type="curve",X=Curve$X,id=Curve$id,time=Curve$time)
  }
  if (is.null(Scalar)==FALSE){
    Scalar <- list(type="scalar",X=Scalar$X,id=Scalar$id)
  }
  if (is.null(Factor)==FALSE){
    Factor <- list(type="factor",X=Factor$X,id=Factor$id)
  }
  
  if (is.null(Shape)==FALSE){
    Shape <- list(type="shape",X=Shape$X,id=Shape$id)
    for (k in 1:dim(Shape$X)[length(dim(Shape$X))]){
      Shape$X[,,,k] <- gpagen(Shape$X[,,,k], print.progress = FALSE)$coords
    }
  }
  
  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs
  
  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }
  #
  
  if (Y$type=="shape"){
    Y$Y <- gpagen(Y$Y,print.progress = FALSE)$coords
  }
  
  # On recupère le nombre de variables au total :
  nvar <- 0
  for (k in Inputs){
    nvar <- nvar + dim(get(k)$X)[length(dim(get(k)$X))]
  }
  
  if (is.null(mtry)==TRUE || mtry> nvar){
    mtry <- floor(nvar/3)*(floor(nvar/3)>=1) + 1*(floor(nvar/3)<1)
  }
  
  if (is.null(Shape)!=TRUE || is.null(Image)!=TRUE) ERT <- TRUE
  
  size <- NULL
  if (Y$type=="shape" || Y$type=="image" || Y$type=="distribution" || Y$type=="covariance" || Y$type=="sphere" || Y$type=="sphere2"){
    size <- c(dim(Y$Y)[1],dim(Y$Y)[2])
  }
  
  if(is.null(ncores)==TRUE){
    ncores <- detectCores()-1
  }
  
  
  print("Building the maximal Frechet trees...")
  
  debut <- Sys.time()
  rf <-  rf_shape_para(Curve=Curve,Scalar=Scalar, Factor=Factor, Shape=Shape, Image=Image,Y=Y, mtry=mtry, ntree=ntree,ERT=ERT,ntry = ntry,timeScale = timeScale,ncores=ncores, aligned.shape = TRUE)
  temps <- Sys.time() - debut
  
  if (Y$type=="shape" || Y$type=="image" || Y$type=="distribution"|| Y$type=="covariance" || Y$type=="sphere" || Y$type=="sphere2"){
    rf <- list(type=Y$type, rf=rf, size = dim(Y$Y) )
  }
  else {
    rf <- list(type=Y$type, rf=rf, levels=levels(Y$Y))
  }
  
  
  print("Forest constucted !")
  xerror <- rep(NA, ntree)
  
  if (OOB==TRUE) {
    print("OOB error computation...")
    
    
    #cl <- parallel::makeCluster(ncores)
    #doParallel::registerDoParallel(cl)
    
    #xerror <- pbsapply(1:ntree, FUN=function(i){OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor,Shape=Shape,Image=Image, Y=Y_KM, timeScale=timeScale,d_out=d_out)},cl=cl)
    
    #parallel::stopCluster(cl)
    xerror = rep(NA,ntree)
    
    if (Y$type=="surv"){
      ykm = NULL
      tkm=NULL
      idkm = NULL
      
      for (i in unique(Y$id)){
        qui = which(Y$id==i)
        donnees = survfit(Surv(Y$time[qui],Y$Y[qui])~1)
        ykm=c(ykm, donnees$surv)
        tkm=c(tkm, donnees$time)
        idkm = c(idkm, rep(i,length(donnees$surv)))
      }
      Y= list(type="surv", Y=ykm, id=idkm, time=tkm)
    }
    
    
    if (Y$type=="surv"){
      for (i in 1:ntree){
        xerror[i] = OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor,Shape=Shape,Image=Image, Y, timeScale=timeScale,d_out=d_out)
      }
      print("Computing the OOB error of the Frechet forest")
      oob.err <- OOB.rfshape(rf,Curve = Curve,Scalar =Scalar,Factor=Factor,Shape=Shape,Image=Image,Y, timeScale=timeScale, d_out=d_out)
    }
    
    else {
      for (i in 1:ntree){
        xerror[i] = OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor,Shape=Shape,Image=Image, Y=Y, timeScale=timeScale,d_out=d_out)
      }
      print("Computing the OOB error of the Frechet forest")
      oob.err <- OOB.rfshape(rf,Curve = Curve,Scalar =Scalar,Factor=Factor,Shape=Shape,Image=Image,Y=Y, timeScale=timeScale, d_out=d_out)
    }
  }
  
  # Ok pour le XERROR
  if (OOB==FALSE && imp == FALSE && Y$type!="surv"){
    frf <- list(rf=rf$rf,type=rf$type,levels=rf$levels, size=size, time=temps)
    class(frf) <- c("FrechForest")
    return(frf)
  }
  
  if (OOB==FALSE && imp == FALSE && Y$type=="surv"){
    frf <- list(rf=rf$rf,type=rf$type,levels=rf$levels, size=size, time=temps)
    class(frf) <- c("FrechForest")
    return(frf)
  } 
  
  if (OOB==TRUE && imp == FALSE && Y$type!="surv"){
    var.ini <- impurity(Y, timeScale=timeScale)
    varex <- 1 - mean(oob.err$err, na.rm = TRUE)/var.ini
    frf <- list(rf=rf$rf,type=rf$type,levels=rf$levels, xerror=xerror,oob.err=oob.err$err,oob.pred= oob.err$oob.pred, varex=varex, size=size, time=temps)
    class(frf) <- c("FrechForest")
    return(frf)
  }
  
  if (OOB==TRUE && imp == FALSE && Y$type=="surv"){
    frf <- list(rf=rf$rf,type=rf$type,levels=rf$levels, xerror=xerror,oob.err=oob.err$err,oob.pred= oob.err$oob.pred, size=size, time=temps)
    class(frf) <- c("FrechForest")
    return(frf)
  }
  
  if (imp == TRUE){
    print("Importance calculation...")
    debut <- Sys.time()
    Curve.perm <- Curve
    Scalar.perm <- Scalar
    Factor.perm <- Factor
    Shape.perm <- Shape
    Image.perm <- Image
    
    Importance.Curve <- NULL
    Importance.Scalar <- NULL
    Importance.Factor <- NULL
    Importance.Shape <- NULL
    Importance.Image <- NULL
    
    #X.perm <- list(type=X$type, X=X$X, id=X$id, time=X$time)
    if (is.element("curve",inputs)==TRUE){
      p=1
      print('Computing the importance on the space of curves')
      Curve.err <- matrix(NA, ntree, dim(Curve$X)[2])
      
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      
      Importance.Curve <- foreach::foreach(p=1:dim(Curve$X)[2],.packages = "kmlShape" ,.combine = "c") %dopar% {
        for (k in 1:ntree){
          BOOT <- rf$rf[,k]$boot
          nboot <- length(unique(Y$id))- length(BOOT)
          
          id_boot_Curve <- NULL
          for (i in 1:length(BOOT)){
            id_boot_Curve <- c(id_boot_Curve, which(Curve$id==BOOT[i]))
          }
          
          # Il faut maintenant faire la permutation :
          
          Curve.perm$X[-id_boot_Curve,p] <- permutation_courbes(Curve$X[-id_boot_Curve,p], Curve$id[-id_boot_Curve])
          
          
          Curve.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve.perm, Scalar = Scalar, Factor=Factor,Shape=Shape, Image=Image, Y, timeScale=timeScale)
          
        }
        Curve.perm$X[,p] <- Curve$X[,p]
        res <- mean(Curve.err[,p]- xerror)
      }
      
      parallel::stopCluster(cl)
    }
    
    
    if (is.element("scalar",inputs)==TRUE){
      p=1
      print('Computing the importance on the space of scalars')
      Scalar.err <- matrix(NA, ntree, dim(Scalar$X)[2])
      
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      
      Importance.Scalar <- foreach::foreach(p=1:dim(Scalar$X)[2],.packages = "kmlShape" ,.combine = "c") %dopar% {
        
        for (k in 1:ntree){
          BOOT <- rf$rf[,k]$boot
          nboot <- length(unique(Y$id))- length(BOOT)
          
          id_boot_Scalar <- NULL
          for (i in 1:length(BOOT)){
            id_boot_Scalar <- c(id_boot_Scalar, which(Scalar$id==BOOT[i]))
          }
          
          
          Scalar.perm$X[-id_boot_Scalar,p] <- sample(Scalar.perm$X[-id_boot_Scalar,p])
          
          Scalar.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar.perm, Factor=Factor,Shape=Shape, Image=Image, Y, timeScale=timeScale)
          
        }
        Scalar.perm$X[,p] <- Scalar$X[,p]
        res <- mean(Scalar.err[,p]- xerror)
      }
      
      parallel::stopCluster(cl)
    }
    
    if (is.element("factor",inputs)==TRUE){
      p=1
      print('Computing the importance on the space of factors')
      Factor.err <- matrix(NA, ntree, dim(Factor$X)[2])
      
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      
      Importance.Factor <- foreach::foreach(p=1:dim(Factor$X)[2],.packages = "kmlShape" ,.combine = "c") %dopar% {
        
        for (k in 1:ntree){
          BOOT <- rf$rf[,k]$boot
          nboot <- length(unique(Y$id))- length(BOOT)
          
          id_boot_Factor <- NULL
          for (i in 1:length(BOOT)){
            id_boot_Factor <- c(id_boot_Factor, which(Factor$id==BOOT[i]))
          }
          
          # Il faut maintenant faire la permutation :
          
          Factor.perm$X[-id_boot_Factor,p] <- sample(Factor.perm$X[-id_boot_Factor,p])
          
          Factor.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar, Factor=Factor.perm ,Shape=Shape, Image=Image, Y, timeScale=timeScale)
          
        }
        ##on remet la variable en place :::
        Factor.perm$X[,p] <- Factor$X[,p]
        res <- mean(Factor.err[,p]- xerror)
      }
      
      parallel::stopCluster(cl)
    }
    
    if (is.element("shape",inputs)==TRUE){
      p=1
      print('Computing the importance on the space of shapes')
      Shape.err <- matrix(NA, ntree, dim(Shape$X)[length(dim(Shape$X))])
      
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      
      Importance.Shape <- foreach::foreach(p=1:dim(Shape$X)[length(dim(Shape$X))],.packages = "kmlShape" ,.combine = "c") %dopar% {
        
        for (k in 1:ntree){
          BOOT <- rf$rf[,k]$boot
          nboot <- length(unique(Y$id))- length(BOOT)
          
          id_boot_Shape <- NULL
          for (i in 1:length(BOOT)){
            id_boot_Shape <- c(id_boot_Shape, which(Shape$id==BOOT[i]))
          }
          
          # Il faut maintenant faire la permutation :
          
          Shape.perm$X[,,-id_boot_Shape,p] <- permutation_shapes(Shape.perm$X[,,-id_boot_Shape, p], Shape.perm$id[-id_boot_Shape])
          
          Shape.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar, Factor=Factor,Shape=Shape.perm, Image=Image, Y, timeScale=timeScale)
          
        }
        ##on remet la variable en place :::
        Shape.perm$X[,,,p] <- Shape$X[,,,p]
        res <- mean(Shape.err[,p]- xerror)
      }
      
      parallel::stopCluster(cl)
    }
    
    if (is.element("image",inputs)==TRUE){
      p=1
      print('Computing the importance on the space of images')
      Image.err <- matrix(NA, ntree, dim(Image$X)[3])
      
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      
      Importance.Image <- foreach::foreach(p=1:dim(Image$X)[3],.packages = "kmlShape" ,.combine = "c") %dopar% {
        
        for (k in 1:ntree){
          BOOT <- rf$rf[,k]$boot
          nboot <- length(unique(Y$id))- length(BOOT)
          
          id_boot_Image <- NULL
          for (i in 1:length(BOOT)){
            id_boot_Image <- c(id_boot_Image, which(Image$id==BOOT[i]))
          }
          
          # Il faut maintenant faire la permutation :
          
          Image.perm$X[-id_boot_Image,,p] <- Image.perm$X[-id_boot_Image,,p][sample(nboot),]
          
          Image.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar, Factor=Factor,Shape=Shape, Image=Image.perm, Y, timeScale=timeScale)
          
        }
        ##on remet la variable en place :::
        Image.perm$X[,,p] <- Image$X[,,p]
        res <- mean(Image.err[,p]- xerror)
      }
      
      parallel::stopCluster(cl)
    }
    
    Importance <- list(Curve=as.vector(Importance.Curve), Scalar=as.vector(Importance.Scalar), Factor=as.vector(Importance.Factor), Shape=as.vector(Importance.Shape), Image=as.vector(Importance.Image))
    
    temps.imp <- Sys.time() - debut
    
    if (OOB==TRUE && Y$type == "surv"){
      frf <- list(rf=rf$rf,type=rf$type,levels=rf$levels,xerror=xerror,oob.err=oob.err$err,oob.pred= oob.err$oob.pred, Importance=Importance, time=temps, size=size)
      class(frf) <- c("FrechForest")
      return(frf)
    }
    if (OOB==TRUE && Y$type != "surv"){
      var.ini <- impurity(Y, timeScale)
      varex <- 1 - mean(oob.err$err)/var.ini
      frf <- list(rf=rf$rf,type=rf$type,levels=rf$levels,xerror=xerror,oob.err=oob.err$err,oob.pred= oob.err$oob.pred, Importance=Importance, varex=varex, time=temps, size=size)
      class(frf) <- c("FrechForest")
      return(frf)
    }
    if (OOB==FALSE && Y$type == "surv"){
      frf <- list(rf=rf$rf,type=rf$type,levels=rf$levels, Importance=Importance, time=temps, size=size)
      class(frf) <- c("FrechForest")
      return(frf)
    }
    if (OOB==FALSE && Y$type != "surv"){
      frf <- list(rf=rf$rf,type=rf$type,levels=rf$levels, Importance=Importance, varex=varex, time=temps, size=size)
      class(frf) <- c("FrechForest")
      return(frf)
    }
    
  }
}
