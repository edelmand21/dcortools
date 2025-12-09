
dcovterms.smp <- function(terms, smp, alg.fast, alg.memsave, alg.standard, p, q, metr.X, metr.Y) {
  
  if (alg.fast) {
    if (p == 1 & q == 1) {
      if (metr.X == "euclidean" & metr.Y == "euclidean") {
          terms.smp <- function(terms, smp) {sampleterms.fast(terms, smp)}
      } else if (metr.X == "euclidean" & metr.Y == "discrete") {
        Y <- as.factor(Y)
        terms.smp <- function(terms, smp) {sampleterms.fast.numdisc(terms, smp)}
      } else if (metr.X == "discrete" & metr.Y == "euclidean") {
        X <- as.factor(X)
        terms.smp <- function(terms, smp) {sampleterms.fast.numdisc(terms, smp)}
      } else if (metr.X == "discrete" & metr.Y == "discrete") {
        X <- as.factor(X)
        Y <- as.factor(Y)
        terms.smp <- function(terms, smp) {sampleterms.fast.discrete(terms, smp)}
      } else {
        stop("metr.X and metr.Y have to be \"euclidean\" or \"discrete\" for fast algorithms")
      }
    } else {
      stop("Dimensions of X and Y must be 1 for fast algorithms.")
    }
  } else if (alg.memsave) {
    if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") &
        metr.Y[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete")) {
      terms.smp <- function(terms, smp) {sampleterms.memsave(terms, smp)}
    } else {
      stop("Memory efficient algorithms cannot be run with user-defined metrics")
    }
  } else if (alg.standard) {
    terms.smp <- function(terms, smp) {sampleterms.standard(terms, smp)}
  }
  return(terms.smp)  
}






sampleterms.standard <- function(terms, smp) {
  aijbij <- matrix_prod_sum_sample(terms$distX, terms$distY, smp)
  Sab <- vector_prod_sum_sample(terms$aidot, terms$bidot, smp)
  return(list("aijbij" = aijbij, "Sab" = Sab))
}


sampleterms.fast.matr <- function(terms, smp) {

 if (terms$ndisc == 0) 
   output <- sampleterms.fast(terms, smp)
 else if (terms$ndisc == 1)
   output <- sampleterms.fast.numdisc(terms, smp)
 else if (terms$ndisc == 2) {
    terms$aidot <- terms$aidotshort
    terms$bidot <- terms$bidotshort
    output <- sampleterms.fast.discrete(terms, smp)
 }
}



sampleterms.fast <- function(terms, smp) {
  IY0.s <- temp <-  1:terms$n
  Sab <- vector_prod_sum_sample(terms$aidot, terms$bidot, smp)
  
  Y.sortX <- terms$Y[smp][terms$IX0]
  XY.sortX <- Y.sortX * terms$X.sort.X
  
  IY.s <- terms$IY[smp][terms$IX0]
  IY0.s[IY.s] <- temp
  
  
  sY.X <- cumsum(terms$X.sort.X[IY0.s])
  sX.Y <- cumsum(Y.sortX)
  sY.XY <- cumsum(XY.sortX[IY0.s]) 
  sX.XY <- cumsum(XY.sortX)
  
  aijbij <- SUMAIJBIJ(IY.s, terms$X.sort.X, Y.sortX, XY.sortX, terms$sX.X, sX.Y, sX.XY, sY.X[IY.s], terms$sY.Y[IY.s], sY.XY[IY.s])
  
  return(list("aijbij" = aijbij, "Sab" = Sab))
}


sampleterms.fast.discrete <- function(terms, smp) {
  
  nXY <- table(terms$X,terms$Y[smp])
  aijbij <- terms$n * (terms$n - 1) - sum(terms$nX * (terms$nX - 1)) - sum(terms$nY * (terms$nY - 1)) + sum(nXY * (nXY - 1))
  Sab <- sum((terms$aidot %*% t(terms$bidot)) * nXY)
  

  return(list("aijbij" = aijbij, "Sab" = Sab))
  
}


sampleterms.fast.numdisc <- function(terms,smp) {
  aijbij <- terms$adotdot
  Y.sortX <- terms$Y.sortX[smp]
  
  
  for (lvl in terms$levY) {
    set0 <- which(Y.sortX == lvl)
    n0 <- length(set0)
    alphaX0 <- 0:(n0-1)
    X0 <- terms$X.sort.X[set0]
    sX.X0 <- cumsum(X0)
    betaX0 <- sX.X0 - X0
    X0dot <- sX.X0[n0]
    aidot0 <- X0dot + (2 * alphaX0 - n0) * X0 - 2 * betaX0
    adotdot0 <- sum(aidot0)
    aijbij <- aijbij - adotdot0
  }
  
  
  bidot <- terms$bidot[smp]

  Sab <- vector_prod_sum(terms$aidot,bidot)
  
  
  return(list("aijbij" = aijbij, "Sab" = Sab))
}


sampleterms.memsave <-  function(terms, smp) {
  

  if (terms$p == 1 & terms$q == 1)
    aijbij <- aijbijmemvec(terms$X, terms$Y[smp], terms$metr.X, terms$metr.Y, terms$prmX, terms$prmY)
    else if (terms$q ==1)
    aijbij <- aijbijmem(terms$X, as.matrix(terms$Y[smp,]), terms$metr.X, terms$metr.Y, terms$prmX, terms$prmY)
    else 
      aijbij <- aijbijmem(terms$X, terms$Y[smp,], terms$metr.X, terms$metr.Y, terms$prmX, terms$prmY)
    
    Sab <- vector_prod_sum_sample(terms$aidot,terms$bidot,smp)
    
    return(list("aijbij" = aijbij, "Sab" = Sab))
}  

