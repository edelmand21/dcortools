prep.fast <- function(X, n, discrete, pairwise) {
  
  cc <- mc <- NULL
  ncc <- n
  
  if (!discrete) {
    IX <- temp <- 1:n
   
    
    if (pairwise) {
      IX0 <- NULL
      cc <- which(stats::complete.cases(X))
      mc <- setdiff(1:n,cc)
      
      ncc <- length(cc)
      IX0cc <- Rfast::Order(X[cc])
      IX[cc][IX0cc] <- 1:ncc
      IX[mc] <- NA
      
      IXcc <- IX[cc]
      Xcc <- X[cc]
      X.sort.X <- Xcc[IX0cc]
    } else {
      IX0 <- Rfast::Order(X)
      X.sort.X <- X[IX0]
      IX[IX0] <- temp
      IXcc <- IX
      Xcc <- X
    }
      
    
    
    sX.X <- cumsum(X.sort.X)
    alphaX <- IXcc - 1
    betaX <- sX.X[IXcc] - X.sort.X[IXcc]
    Xdot <- sX.X[ncc]
    aidot <- Xdot + (2 * alphaX - ncc) * Xcc - 2 * betaX
    adotdot <- sum(aidot)
    
    XXdot <- sum(Xcc^2)
    aijaij <- 2 * ncc * XXdot - 2* Xdot^2
    
    return(list("cc" =cc, "mc" = mc, "IX"= IX, "IX0" = IX0, "X.sort.X" = X.sort.X, "aidot" = aidot, "adotdot" = adotdot, "aijaij" = aijaij, "ncc" = ncc, "sX.X" = sX.X, "X" = X))
  } else {
    Xcc <- X
    if (pairwise) {
      cc <- which(stats::complete.cases(X))
      n <- length(cc)
      Xcc <- X[cc]
    }
    nX <- table(Xcc)
    aijaij <- n * (n - 1) - sum(nX * (nX - 1)) 
    aidots <- as.numeric(n - nX)
    adotdot <- sum(aidots * nX)
    aidot <- as.numeric(n - nX[Xcc])
    
    return(list("cc" =cc, "aidotshort" = aidots, "aidot" = aidot, "adotdot" = adotdot, "aijaij" = aijaij, "X" = as.factor(X), "nX" = nX, "ncc" = n))
  }
}





prep.standard <- function(X, n, p,  metr.X, pairwise) {
  
  distX <- distmat(X, metr.X, p)
  
  if (!pairwise) {
    aidot <- Rfast::colsums(distX)
    adotdot <- sum(aidot)
    aijaij <- matrix_prod_sum(distX, distX)
    
    return(list("distX" = distX, "aidot" = aidot, "adotdot" = adotdot, "aijaij" = aijaij, "ncc" = n))
    
  } else {
    cc <- which(stats::complete.cases(X))
    ncc <- length(cc)
    aidot <- colsums_subset(distX,cc)
    adotdot <- sum(aidot)
    aijaij <- matrix_prod_sum_subset(distX, distX, cc)
    return(list("distX" = distX, "aidot" = aidot, "adotdot" = adotdot, "aijaij" = aijaij,  "ncc" = ncc, "cc" = cc))
  }
  
}

prep.memsave <- function(X, n, p,  metr.X, pairwise) {
  if (pairwise) {
    cc <- which(stats::complete.cases(X))
    Xcc <- X[cc]
    n <- length(cc)
  } else {
    Xcc <- X
  }
  
  output <- dvarterms.memsave(Xcc, metr.X, p)
  
  output$ncc <- n
  output$p <- p
  output$X <- X
  
  if (pairwise) {
    output$cc <- cc
  }
  
  
  return(output)
  
}
  