
dcovterms <- function(X,Y,n, calc.dcor, doperm, dobb3, alg.fast, alg.memsave, alg.standard, p, q, metr.X, metr.Y, type.X, type.Y) {

    if (alg.fast) {
    if (p == 1 & q == 1) {
      if (metr.X == "euclidean" & metr.Y == "euclidean") {
        terms <- dcovterms.fast(X, Y, n, calc.dcor = calc.dcor, calc.perm = doperm)
      } else if (metr.X == "euclidean" & metr.Y == "discrete") {
        Y <- as.factor(Y)
        terms <- dcovterms.fast.numdisc(X, Y, n, calc.dcor = calc.dcor, calc.perm = doperm)
      } else if (metr.X == "discrete" & metr.Y == "euclidean") {
        X <- as.factor(X)
        terms <- dcovterms.fast.numdisc(Y, X, n, calc.dcor = calc.dcor, calc.perm = doperm)
      } else if (metr.X == "discrete" & metr.Y == "discrete") {
        X <- as.factor(X)
        Y <- as.factor(Y)
        terms <- dcovterms.fast.discrete(X ,Y, n, calc.dcor = calc.dcor, calc.perm = doperm)
      } else {
        stop("metr.X and metr.Y have to be \"euclidean\" or \"discrete\" for fast algorithms")
      }
    } else {
      stop("Dimensions of X and Y must be 1 for fast algorithms.")
    }
  } else if (alg.memsave) {
    if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") &
        metr.Y[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete")) {
      terms <- dcovterms.memsave(X, Y, metr.X, metr.Y, p, q, calc.dcor = calc.dcor, calc.perm = doperm)
    } else {
      stop("Memory efficient algorithms cannot be run with user-defined metrics")
    }
  } else if (alg.standard) {
    terms <- dcovterms.standard(X, Y, type.X, type.Y, metr.X, metr.Y, p, q, calc.dcor = calc.dcor, calc.bb3 = dobb3 | doperm)
}
 return(terms)  
}






dcovterms.standard <- function(X, Y = NULL, type.X, type.Y = NULL, metr.X, metr.Y = NULL, p, q = NULL, calc.dvar = FALSE, calc.dcor = FALSE, calc.bb3 = FALSE) {
  ## if distance matrix is given
  if (type.X == "distance") {
    distX <- X
  } else {
    distX <- distmat(X, metr.X, p)
  }
  
  aidot <- Rfast::colsums(distX)
  adotdot <- sum(aidot)
  
  if (calc.dvar) {
    aijaij <- matrix_prod_sum(distX, distX)
    Saa <- vector_prod_sum(aidot, aidot)
    Taa <- adotdot * adotdot
    return(list("aijaij" = aijaij, "Saa" = Saa, "Taa" = Taa))
  }
  
  if (type.Y == "distance") {
    distY <- Y
  } else {
    distY <- distmat(Y, metr.Y, q)
  }
  

  bidot <- Rfast::colsums(distY)
  bdotdot <- sum(bidot)
  
  aijbij <- matrix_prod_sum(distX, distY)
  Sab <- vector_prod_sum(aidot, bidot)
  
  if (calc.dcor) {
    aijaij <- matrix_prod_sum(distX, distX)
    bijbij <- matrix_prod_sum(distY, distY)
    Saa <- vector_prod_sum(aidot, aidot)
    Sbb <- vector_prod_sum(bidot, bidot)
    if (calc.bb3)
      return(list("aijbij" = aijbij, "Sab" = Sab, "aijaij" = aijaij, "Saa" = Saa, "adotdot" = adotdot, "bijbij" = bijbij, "Sbb" = Sbb, "bdotdot" = bdotdot, "aidot" = aidot, "bidot" = bidot, "distX" = distX, "distY" = distY))
      else
      return(list("aijbij" = aijbij, "Sab" = Sab, "aijaij" = aijaij, "Saa" = Saa, "adotdot" = adotdot, "bijbij" = bijbij, "Sbb" = Sbb, "bdotdot" = bdotdot))
  }
  

  
  Tab <- adotdot * bdotdot
  
  return(list("aijbij" = aijbij, "Sab" = Sab, "Tab" = Tab))
}



dcovterms.fast <- function(X, Y = NULL, n, calc.dvar = FALSE, calc.dcor = FALSE, calc.perm = FALSE) {
  X <- as.numeric(X)
  Y <- as.numeric(Y)
  
  temp <- IX <- IY <- IY.s <- 1:n
  #
  
  IX0 <- Rfast::Order(X)
  X.sort.X <- X[IX0]
  IX[IX0] <- temp
  
  sX.X <- cumsum(X.sort.X)
  alphaX <-  0:(n - 1)
  betaX <- sX.X - X.sort.X
  Xdot <- sX.X[n]
  aidot <- Xdot + (2 * alphaX - n) * X.sort.X - 2 * betaX
  adotdot <- sum(aidot)
  
  if (calc.dvar) {
    XXdot <- sum(X ^ 2)
    aijaij <- 2 * n * XXdot - 2* Xdot^2
    Saa <- vector_prod_sum(aidot, aidot)
    Taa <- adotdot * adotdot
    return(list("aijaij" = aijaij, "Saa" = Saa, "Taa" = Taa))
  }
  
  IY0 <- Rfast::Order(Y)
  Y.sort.Y <- Y[IY0]
  IY[IY0] <- temp
  alphaY <- 0:(n - 1)
  
  sY.Y <- cumsum(Y.sort.Y)
  
  betaY <- sY.Y - Y.sort.Y
  
  Ydot <- sY.Y[n]
  
  bidot <- Ydot + (2 * alphaY - n) * Y.sort.Y - 2 * betaY
  bdotdot <- sum(bidot)
  
  Y.sortX <- Y[IX0]
  XY.sortX <- Y.sortX * X.sort.X
  IY0.s <- IX[IY0]
  IY.s[IY0.s] <- temp
  
  
  sY.X <- cumsum(X.sort.X[IY0.s])
  sX.Y <- cumsum(Y.sortX)
  sY.XY <- cumsum(XY.sortX[IY0.s]) 
  sX.XY <- cumsum(XY.sortX)
  
 
  aijbij <- SUMAIJBIJ(IY.s, X.sort.X, Y.sortX, XY.sortX, sX.X, sX.Y, sX.XY, sY.X[IY.s], sY.Y[IY.s], sY.XY[IY.s])
 
  Sab <- vector_prod_sum(aidot, bidot[IY.s])
  
  if (calc.dcor) {
    XXdot <- sum(X^2)
    YYdot <- sum(Y^2)
    Saa <- vector_prod_sum(aidot, aidot)
    aijaij <- 2 * n * XXdot - 2* Xdot^2
    Sbb <- vector_prod_sum(bidot, bidot)
    bijbij <- 2 * n * YYdot - 2* Ydot^2
    if (calc.perm) {
      aidot <- aidot[IX]
      bidot <- bidot[IY]
      return(list("aijbij" = aijbij, "Sab" = Sab, "aijaij" = aijaij, "Saa" = Saa, "adotdot" = adotdot, "bijbij" = bijbij, "Sbb" = Sbb, "bdotdot" = bdotdot, 
                  "aidot" = aidot, "bidot" = bidot, "X.sort.X" = X.sort.X, "Y" = Y, "IX0" = IX0, "sX.X" = sX.X, "sY.Y" = sY.Y, "IY" = IY, "n" = n, "ndisc"=0))
    }
      else
      return(list("aijbij" = aijbij, "Sab" = Sab, "aijaij" = aijaij, "Saa" = Saa, "adotdot" = adotdot, "bijbij" = bijbij, "Sbb" = Sbb, "bdotdot" = bdotdot))
  }
  
  
  Tab <- adotdot * bdotdot
  
  return(list("aijbij" = aijbij, "Sab" = Sab, "Tab" = Tab))
}


dcovterms.fast.discrete <- function(X, Y = NULL, n, calc.dvar = FALSE, calc.dcor = FALSE, calc.perm = FALSE) {
  
  if (calc.dvar) {
    labX <- levels(X)
    fX <- length(labX)
    nX <- table(X)
    aijaij <- n * (n - 1) - sum(nX * (nX - 1)) 
    aidot <- as.numeric(n - nX)
    adotdot <- sum(aidot * nX)
    Saa <- sum(aidot * aidot * nX)
    Taa <- adotdot*adotdot
    return(list("aijaij" = aijaij, "Saa" = Saa, "Taa" = Taa))
  } 
  
  nXY <- Rfast::Table(X, Y, names = FALSE)
  nX <- Rfast::rowsums(nXY)
  nY <- Rfast::colsums(nXY)
  
  
  aidot <- as.numeric(n - nX)
  adotdot <- sum(aidot * nX)
  bidot <- as.numeric(n - nY)
  bdotdot <- sum(bidot * nY)
  
  aijbij <- n * n - sum(nX * nX) - sum(nY * nY) + sum(nXY * nXY)
  Sab <- sum((aidot %*% t(bidot)) * nXY)
  
  
  if (calc.dcor) {
    aijaij <- n * (n - 1) - sum(nX * (nX - 1)) 
    bijbij <- n * (n - 1) - sum(nY * (nY - 1)) 
    Saa <- sum(aidot * aidot * nX)
    Sbb <- sum(bidot * bidot * nY)
    if (calc.perm) {
      return(list("aijbij" = aijbij, "Sab" = Sab, "aijaij" = aijaij, "Saa" = Saa, "adotdot" = adotdot, "bijbij" = bijbij, "Sbb" = Sbb, "bdotdot" = bdotdot,
                  "X" = X, "Y" = Y, "nX" = nX, "nY" = nY, "aidot" = aidot, "bidot" = bidot, "n" = n,"ndisc"=2)) 
    } else {
      return(list("aijbij" = aijbij, "Sab" = Sab, "aijaij" = aijaij, "Saa" = Saa, "adotdot" = adotdot, "bijbij" = bijbij, "Sbb" = Sbb, "bdotdot" = bdotdot))
    }
  } 
  
  Tab <- adotdot*bdotdot
  
  return(list("aijbij" = aijbij, "Sab" = Sab, "Tab" = Tab))
  
}


dcovterms.fast.numdisc <- function(X, Y, n, calc.dcor =FALSE, calc.perm = FALSE) {
  IX0 <- Rfast::Order(X)
  X.sort.X <- X[IX0]

  
  sX.X <- cumsum(X.sort.X)
  alphaX <- 0:(n - 1)
  betaX <- sX.X - X.sort.X
  Xdot <- sX.X[n]
  Y.sortX <- Y[IX0]
  
  aidot <- Xdot + (2 * alphaX - n) * X.sort.X - 2 * betaX
  adotdot <- sum(aidot)
  levY <- levels(Y)
  nY <- Rfast::Table(Y, names = FALSE)
  aijbij <- adotdot
  
  
  for (lvl in levY) {
    set0 <- which(Y.sortX == lvl)
    n0 <- length(set0)
    alphaX0 <- 0:(n0 - 1)
    X0 <- X.sort.X[set0]
    sX.X0 <- cumsum(X0)
    betaX0 <- sX.X0 - X0
    X0dot <- sX.X0[n0]
    aidot0 <- X0dot + (2 * alphaX0 - n0) * X0 - 2 * betaX0
    adotdot0 <- sum(aidot0)
    aijbij <- aijbij - adotdot0
  }
  
  bidot <- as.numeric(n - nY[Y.sortX])
  bdotdot <- sum(bidot)
  
  Sab <- vector_prod_sum(aidot, bidot)
  
  
  if (calc.dcor) {
    XXdot <- sum(X ^ 2)
    aijaij <- 2 * n * XXdot - 2 * Xdot ^ 2
    bijbij <- n * (n - 1) - sum(nY * (nY - 1))
    Saa <- sum(aidot * aidot)
    Sbb <- sum(bidot * bidot)
    if (calc.perm)
      return(
        list(
          "aijbij" = aijbij,
          "Sab" = Sab,
          "aijaij" = aijaij,
          "Saa" = Saa,
          "adotdot" = adotdot,
          "bijbij" = bijbij,
          "Sbb" = Sbb,
          "bdotdot" = bdotdot,
          "aidot" = aidot,
          "bidot" = bidot,
          "X.sort.X" = X.sort.X,
          "Y.sortX" = Y.sortX,
          "levY" = levY,
          "n" = n,
          "ndisc"=1
        )
      )
    else
      return(
        list(
          "aijbij" = aijbij,
          "Sab" = Sab,
          "aijaij" = aijaij,
          "Saa" = Saa,
          "adotdot" = adotdot,
          "bijbij" = bijbij,
          "Sbb" = Sbb,
          "bdotdot" = bdotdot
        )
      )
  } 
  
  
  Tab <- adotdot * bdotdot
  
  return(list("aijbij" = aijbij, "Sab" = Sab, "Tab" = Tab))
}


dcovterms.memsave <-  function(X, Y, metr.X, metr.Y, p, q, calc.dcor = FALSE, calc.perm = FALSE) {
  
  if (length(metr.X) == 1) {
    prmX <- 0
  } else {
    prmX <- as.numeric(metr.X[2])
    metr.X <- as.character(metr.X[1])
    if (is.na(prmX)) 
      stop("Parameters of standard metrics must be numeric")
  }
  
  if (length(metr.Y) == 1) {
    prmY <- 0
  } else {
    prmY <- as.numeric(metr.Y[2]) 
    metr.Y <- metr.Y[1]
    if (is.na(prmY)) 
      stop("Parameters of standard metrics must be numeric")
  }
  
  
  
  if (p == 1 & q == 1) {
    if (metr.X == "discrete") {
      X <- as.numeric(X)
    }
    if (metr.Y == "discrete") {
      Y <- as.numeric(Y)
    }
    terms <- dcovtermsmemvec(X, Y, metr.X, metr.Y, prmX, prmY, calcdcor = calc.dcor, calcperm = calc.perm)
  } else if (p > 1 & q > 1) {
    terms <- dcovtermsmem(X, Y, metr.X, metr.Y, prmX, prmY, calcdcor = calc.dcor, calcperm = calc.perm)
  } else if (p > 1) {
    if (metr.Y == "discrete") {
      Y <- as.numeric(Y)
    }
    Y <- as.matrix(Y)
    terms <- dcovtermsmem(X, as.matrix(Y), metr.X, metr.Y, prmX, prmY, calcdcor = calc.dcor, calcperm = calc.perm)
  } else {
    if (metr.X == "discrete") {
      X <- as.numeric(X)
    }
    X <- as.matrix(X)
    terms <- dcovtermsmem(X, Y, metr.X, metr.Y, prmX, prmY, calcdcor = calc.dcor, calcperm = calc.perm)
  }
  
  if (calc.perm) {
    terms$prmX <- prmX
    terms$prmY <- prmY
    terms$metr.X <- metr.X
    terms$metr.Y <- metr.Y
    terms$p <- p
    terms$q <- q
    terms$X <- X
    terms$Y <- Y
  }
  
  return(terms)
}  



dvarterms.memsave <-  function(X, metr.X, p) {
  
  if (length(metr.X) == 1) {
    prmX <- 0
  } else {
    prmX <- as.numeric(metr.X[2])
    metr.X <- as.character(metr.X[1])
    if (is.na(prmX)) 
      stop("Parameters of standard metrics must be numeric")
  }
  

  if (p == 1) {
    if (metr.X == "discrete") {
      X <- as.numeric(X)
    }
    terms <- dvartermsmemvec(X, metr.X, prmX)
  } else {
    terms <- dvartermsmem(X, metr.X, prmX)
  }
  
  terms$Saa <- vector_prod_sum(terms$aidot,terms$aidot)
  terms$Taa <- terms$adotdot * terms$adotdot
  
  return(terms)
}  


aijbij.memsave <-  function(X, Y, metr.X, metr.Y, p, q) {
  
  if (length(metr.X) == 1) {
    prmX <- 0
  } else {
    prmX <- as.numeric(metr.X[2])
    metr.X <- as.character(metr.X[1])
    if (is.na(prmX)) 
      stop("Parameters of standard metrics must be numeric")
  }
  
  if (length(metr.Y) == 1) {
    prmY <- 0
  } else {
    prmY <- as.numeric(metr.Y[2]) 
    metr.Y <- metr.Y[1]
    if (is.na(prmY)) 
      stop("Parameters of standard metrics must be numeric")
  }
  

  if (p == 1 & q == 1) {
    if (metr.X == "discrete") {
      X <- as.numeric(X)
    }
    if (metr.Y == "discrete") {
      Y <- as.numeric(Y)
    }
    terms <- aijbijmemvec(X, Y, metr.X, metr.Y, prmX, prmY)
  } else if (p > 1 & q > 1) {
    terms <- aijbijmem(X, Y, metr.X, metr.Y, prmX, prmY)
  } else if (p > 1) {
    if (metr.Y == "discrete") {
      Y <- as.numeric(Y)
    }
    Y <- as.matrix(Y)
    terms <- aijbijmem(X, Y, metr.X, metr.Y, prmX, prmY)
  } else {
    if (metr.X == "discrete") {
      X <- as.numeric(X)
    }
    X <- as.matrix(X)
    terms <- aijbijmem(X, Y, metr.X, metr.Y, prmX, prmY)
  }
  
  
  return(terms)
}  
