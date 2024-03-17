#########################revised versio of IntrinsicMean from #########################
IntrinsicMean <- function(data, weights = rep(1, nrow(data)), thres = 1e-5){         # 'data': longitude/latitude matrix or dataframe with two column.
  if (sum(weights^2)^(1/2) < 1e-15) stop("weight should not be the zero vector.")     # 'weights': n-dimensional vector.
  t <- 0
  mu <- data[1, , drop = F]                                                          # Initialize mean as a point.
  delta.mu <- c(1, 0)
  while (sum((delta.mu)^2)^(1/2) > thres){
    weights.normal <- weights/sum(weights)                                           # Normalize of 'weights' so that its components add up to 1.
    summation <- c(0, 0)
    for (m in 1:length(weights)){
      rot <- Rotate(mu, data[m, ])
      summation <- summation + weights.normal[m] * Logmap(rot)
    }
    delta.mu <- summation
    exp <- Expmap(delta.mu)
    mu.Euc <- Rotate.inv(mu, Trans.sph(exp))
    mu <- Trans.sph(mu.Euc)
    t <- t+1
    if (t>100) {break}
  }
  if (t<=100){
    return(as.numeric(mu))
  }else{
    return(c(5,5))
  }
}