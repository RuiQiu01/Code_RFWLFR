##################log_cholesky度量##########################
distLog_cholesky <- function(A, B) {
  M.a <- chol(A)
  D.a <- diag(M.a)
  L.a <- M.a - diag(D.a)
  M.b <- chol(B)
  D.b <- diag(M.b)
  L.b <- M.b - diag(D.b) 
  sqrt(sum((L.a-L.b)^2)+sum((log(D.a)-log(D.b))^2)) 
}

##################log_choleskey度量下weighted freshet mean################
estLog_cholesky <- function(S, weights = 1) {
  M <- dim(S)[3]
  if (length(weights) == 1) {
    weights <- rep(1, times = M)
  }
  sum <- S[, , 1] * 0
  for (j in 1:M) {
    C <- chol(S[, , j])
    D <- diag(C)
    L <- C - diag(D)
    sum <- sum + t(diag(log(D))+L) * weights[j] / sum(weights)
  }
  cc <- sum
  D <- diag(cc)
  L <- cc - diag(D)
  cc <- diag(exp(D))+L
  cc %*% t(cc)
}

estcov <-
  function (S,
            method = "Riemannian",
            weights = 1,
            alpha = 1 / 2,
            MDSk = 2)
  {
    out <- list(
      mean = 0,
      sd = 0,
      pco = 0,
      eig = 0,
      dist = 0
    )
    M <- dim(S)[3]
    if (length(weights) == 1) {
      weights <- rep(1, times = M)
    }
    if (method == "Procrustes") {
      dd <- estSS(S, weights)
    }
    if (method == "ProcrustesShape") {
      dd <- estShape(S, weights)
    }
    if (method == "Riemannian") {
      dd <- estLogRiem2(S, weights)
    }
    if (method == "Cholesky") {
      dd <- estCholesky(S, weights)
    }
    if (method == "Log_cholesky") {
      dd <- estLog_cholesky(S, weights)
    }
    if (method == "Power") {
      dd <- estPowerEuclid(S, weights, alpha)
    }
    if (method == "Euclidean") {
      dd <- estEuclid(S, weights)
    }
    if (method == "LogEuclidean") {
      dd <- estLogEuclid(S, weights)
    }
    if (method == "RiemannianLe") {
      dd <- estRiemLe(S, weights)
    }
    out$mean <- dd
    out
  }


distcov <- function(S1,
                    S2 ,
                    method = "Riemannian",
                    alpha = 1 / 2) {
   if (method == "Procrustes") {
      dd <- distProcrustesSizeShape(S1, S2)
   }
   if (method == "ProcrustesShape") {
      dd <- distProcrustesFull(S1, S2)
   }
   if (method == "Riemannian") {
      dd <- distRiemPennec(S1, S2)
   }
   if (method == "Cholesky") {
      dd <- distCholesky(S1, S2)
   }
   if (method == "Log_cholesky") {
      dd <- distLog_cholesky(S1, S2)
   }
   if (method == "Power") {
      dd <- distPowerEuclidean(S1, S2, alpha)
   }
   if (method == "Euclidean") {
      dd <- distEuclidean(S1, S2)
   }
   if (method == "LogEuclidean") {
      dd <- distLogEuclidean(S1, S2)
   }
   if (method == "RiemannianLe") {
      dd <- distRiemannianLe(S1, S2)
   }
   dd
}
